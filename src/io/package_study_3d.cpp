#include "nova/io/package_study_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace nova::io {
namespace {
using physics::TimeSignal3D;
struct Record {
    std::size_t line{};
    std::vector<std::string> words;
    void count(std::size_t n) const {
        if (words.size()!=n) throw std::invalid_argument("Wrong number of arguments for "+words.front());
    }
};
template<typename T> T number(const std::string& word) {
    T val{};
    const auto result=std::from_chars(word.data(),word.data()+word.size(),val);
    if (result.ec!=std::errc{} || result.ptr!=word.data()+word.size() || !std::isfinite(val))
        throw std::invalid_argument("Invalid number "+word);
    return val;
}
std::string contents(const std::filesystem::path& path) {
    std::ifstream input(path,std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open "+path.string());
    std::ostringstream out; out<<input.rdbuf();
    if (input.bad() || !out) throw std::runtime_error("Cannot read "+path.string());
    return out.str();
}
// Framed FNV-1a detects accidental edits including mesh tags. Not authentication.
std::string identity(const std::string& model,const std::string& mesh) {
    std::uint64_t hash=14695981039346656037ULL;
    const auto mix=[&](unsigned char c){hash=(hash^c)*1099511628211ULL;};
    for (const auto* bytes:{&model,&mesh}) {
        const auto length=static_cast<std::uint64_t>(bytes->size());
        for (unsigned i=0;i<8;++i) mix(static_cast<unsigned char>(length>>(8*i)));
        for (unsigned char c:*bytes) mix(c);
    }
    std::ostringstream out; out<<std::hex<<hash; return out.str();
}
bool identifier(const std::string& name) {
    return !name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")==std::string::npos;
}
PackageProbe3D locate(const fem::TetrahedralMesh3D& mesh,const std::string& name,fem::Point3D p) {
    for (std::size_t e=0;e<mesh.elements.size();++e) {
        const auto& a=mesh.nodes[mesh.elements[e][0]];
        const auto g=fem::evaluate_tetrahedron_p1_physical_gradients_3d(mesh,e);
        std::array<double,4> weights{};
        for (unsigned j=0;j<4;++j) {
            const auto& v=g.gradients[j];
            weights[j]=(j==0 ? 1. : 0.)+v[0]*(p.x-a.x)+v[1]*(p.y-a.y)+v[2]*(p.z-a.z);
        }
        if (std::all_of(weights.begin(),weights.end(),[](double w){return w>=-1e-10 && w<=1+1e-10;}))
            return {name,e,weights};
    }
    throw std::invalid_argument("Probe lies outside the mesh: "+name);
}
}
PackageStudy3D read_package_study_3d(const std::filesystem::path& path) {
    const auto source=contents(path);
    std::istringstream lines(source);
    std::vector<Record> rec;
    std::string line;
    for (std::size_t no=1;std::getline(lines,line);++no) {
        std::istringstream words(line.substr(0,line.find('#')));
        Record r{no,{}};
        for (std::string word;words>>word;) r.words.push_back(word);
        if (!r.words.empty()) rec.push_back(std::move(r));
    }
    if (rec.size()<2 || rec[0].words!=std::vector<std::string>{"NOVA_PACKAGE","1"} ||
        rec[1].words.size()!=2 || rec[1].words[0]!="mesh")
        throw std::invalid_argument("Expected NOVA_PACKAGE 1 followed by mesh <Gmsh file>.");
    const auto mesh_data=contents(path.parent_path()/rec[1].words[1]);
    std::istringstream mesh_input(mesh_data);
    auto imported=read_gmsh_tetrahedra(mesh_input);
    PackageStudy3D study;
    auto& m=study.model; auto& o=study.options;
    study.model_id=identity(source,mesh_data);
    study.material_tags=imported.volume_physical_tags;
    const std::set<int> volumes(study.material_tags.begin(),study.material_tags.end());
    std::map<int,std::vector<fem::BoundaryFace3D>> boundary;
    for (std::size_t i=0;i<imported.exterior_faces.size();++i)
        boundary[imported.exterior_physical_tags[i]].push_back(imported.exterior_faces[i]);
    std::map<int,std::size_t> materials;
    std::map<int,TimeSignal3D> sources;
    std::map<std::string,TimeSignal3D> curves;
    std::set<std::string> singles,probe_names;
    std::map<std::size_t,std::pair<std::string,int>> voltage, temperature;
    std::set<std::size_t> clamped;
    std::set<std::pair<std::string,int>> face_settings;
    const auto get_curve=[&](const std::string& name) {
        const auto found=curves.find(name);
        if (found==curves.end()) throw std::invalid_argument("Undefined curve "+name);
        return found->second;
    };
    for (std::size_t index=2;index<rec.size();++index) {
        const auto& r=rec[index]; const auto& w=r.words; const auto& key=w[0];
        const auto x=[&](std::size_t i){return number<double>(w.at(i));};
        const auto single=[&](){if (!singles.insert(key).second) throw std::invalid_argument("Repeated "+key);};
        try {
            if (key=="material" || key=="source") {
                r.count(key=="material" ? 13 : 3);
                const auto tag=physical_tag(w[1],3,imported);
                if (!volumes.contains(tag)) throw std::invalid_argument("Unknown material volume "+w[1]);
                if (key=="source") {
                    if (!sources.emplace(tag,get_curve(w[2])).second) throw std::invalid_argument("Duplicate source.");
                } else {
                    if (!materials.emplace(tag,m.materials.size()).second) throw std::invalid_argument("Duplicate material.");
                    m.materials.push_back({x(2),{x(3),x(4)},{x(5),x(6)},{x(7),x(8)},{x(9),x(10)},x(11),x(12)});
                }
            } else if (key=="tensor") {
                r.count(9);
                const auto tag=physical_tag(w[2],3,imported);
                if (!materials.contains(tag)) throw std::invalid_argument("Define material before its tensor.");
                if (w[1]!="thermal" && w[1]!="electrical") throw std::invalid_argument("Unknown tensor field.");
                auto& mat=m.materials[materials.at(tag)];
                auto& shape=w[1]=="thermal" ? mat.thermal_shape : mat.electrical_shape;
                if (shape) throw std::invalid_argument("Repeated material tensor.");
                shape=fem::SymmetricTensor3D{x(3),x(4),x(5),x(6),x(7),x(8)};
                fem::validate_conductivity_tensor_3d(*shape);
            } else if (key=="curve") {
                if (w.size()<4 || w.size()%2!=0 || !identifier(w[1])) throw std::invalid_argument("Expected named time/value pairs.");
                std::vector<std::pair<double,double>> points;
                for (std::size_t j=2;j<w.size();j+=2) {
                    if (!points.empty() && x(j)<=points.back().first) throw std::invalid_argument("Curve times must increase.");
                    points.emplace_back(x(j),x(j+1));
                    if (w.size()>4) o.events.push_back(x(j));
                }
                auto value=[p=std::move(points)](double t) {
                    const auto b=std::upper_bound(p.begin(),p.end(),t,[](double a,const auto& v){return a<v.first;});
                    if (b==p.begin()) return b->second;
                    const auto a=std::prev(b);
                    if (b==p.end()) return a->second;
                    return std::lerp(a->second,b->second,(t-a->first)/(b->first-a->first));
                };
                if (!curves.emplace(w[1],std::move(value)).second) throw std::invalid_argument("Duplicate curve.");
            } else if (key=="voltage" || key=="temperature" || key=="clamp" || key=="cool") {
                if (key=="cool" && w.size()==6) r.count(6);
                else r.count(key=="cool" ? 5 : (key=="clamp" ? 2 : 3));
                const auto tag=physical_tag(w[1],2,imported);
                if (!boundary.contains(tag) || !face_settings.emplace(key,tag).second)
                    throw std::invalid_argument("Unknown or repeated surface "+w[1]);
                for (const auto& face:boundary.at(tag)) {
                    if (key=="cool") {
                        m.cooling.push_back({face,x(2),x(3),get_curve(w[4]),w.size()==6 ? get_curve(w[5]) : TimeSignal3D{}});
                    } else for (auto i:face) {
                        if (key=="clamp") clamped.insert(i);
                        else {
                            static_cast<void>(get_curve(w[2]));
                            auto& constraints=key=="voltage" ? voltage : temperature;
                            const auto [it,inserted]=constraints.emplace(i,std::pair{w[2],tag});
                            if (!inserted && (it->second.first!=w[2] || (key=="voltage" && it->second.second!=tag)))
                                throw std::invalid_argument("Conflicting constraints or touching terminals at a node.");
                        }
                    }
                }
            } else if (key=="initial") { r.count(2); single(); study.initial_kelvin=x(1); }
            else if (key=="study") {
                r.count(6); single();
                if (w[1]=="be") o.scheme=time::ImplicitScheme::BackwardEuler;
                else if (w[1]=="cn") o.scheme=time::ImplicitScheme::CrankNicolson;
                else if (w[1]=="bdf2") o.scheme=time::ImplicitScheme::BDF2;
                else throw std::invalid_argument("Unknown integration scheme.");
                study.end_time=x(2); o.initial_dt=x(3); o.minimum_dt=x(4); o.maximum_dt=x(5);
            } else if (key=="solver") {
                r.count(5); single();
                if (w[1]=="monolithic") o.strategy=physics::TransportStrategy3D::MonolithicNewton;
                else if (w[1]=="segregated") o.strategy=physics::TransportStrategy3D::Segregated;
                else throw std::invalid_argument("Unknown nonlinear strategy.");
                o.newton.residual_tolerance=x(2); o.relative_residual=x(3); o.newton.max_iterations=number<std::size_t>(w[4]);
            } else if (key=="linear") {
                r.count(4); single();
                if (w[1]=="gmres") o.newton.linear_solver=nonlinear::NewtonLinearSolver::GMRES;
                else if (w[1]=="bicgstab") o.newton.linear_solver=nonlinear::NewtonLinearSolver::BiCGSTAB;
                else throw std::invalid_argument("Unknown transport linear solver.");
                linalg::PreconditionerType pre{};
                if (w[2]=="ilu0") pre=linalg::PreconditionerType::IncompleteLU0;
                else if (w[2]=="jacobi") pre=linalg::PreconditionerType::Jacobi;
                else if (w[2]!="none") throw std::invalid_argument("Unknown transport preconditioner.");
                const auto limit=number<std::size_t>(w[3]);
                o.newton.gmres.solver={limit,1e-11,1e-15,pre}; o.newton.bicgstab=o.newton.gmres.solver;
            } else if (key=="adapt") { r.count(3); single(); o.error_kelvin=x(1); o.relative_error=x(2); }
            else if (key=="scales") { r.count(3); single(); o.temperature_scale=x(1); o.voltage_scale=x(2); }
            else if (key=="relaxation") { r.count(2); single(); o.relaxation=x(1); }
            else if (key=="frames") { r.count(2); single(); study.frame_stride=number<std::size_t>(w[1]); }
            else if (key=="probe") {
                r.count(5);
                if (!identifier(w[1]) || !probe_names.insert(w[1]).second) throw std::invalid_argument("Invalid probe name.");
                study.probes.push_back(locate(imported.mesh,w[1],{x(2),x(3),x(4)}));
            } else throw std::invalid_argument("Unknown model directive "+key);
        } catch (const std::exception& e) {
            throw std::invalid_argument(path.string()+" line "+std::to_string(r.line)+": "+e.what());
        }
    }
    if (!singles.contains("study") || !singles.contains("initial") || study.initial_kelvin<=0 || study.end_time<=0)
        throw std::invalid_argument("Model requires positive initial Kelvin temperature and end time.");
    for (auto tag:study.material_tags) {
        if (!materials.contains(tag)) throw std::invalid_argument("Missing material for volume "+std::to_string(tag));
        m.cell_material.push_back(materials.at(tag));
    }
    m.heat_sources.resize(m.materials.size());
    for (const auto& [tag,curve]:sources) m.heat_sources[materials.at(tag)]=curve;
    for (const auto& [node,bc]:voltage) {
        m.voltage.push_back({node,get_curve(bc.first)}); study.terminal_tags.push_back(bc.second);
    }
    for (const auto& [node,bc]:temperature) m.temperature.push_back({node,get_curve(bc.first)});
    for (auto node:clamped) m.clamps.push_back({node,{}});
    std::sort(o.events.begin(),o.events.end()); o.events.erase(std::unique(o.events.begin(),o.events.end()),o.events.end());
    m.mesh=std::move(imported.mesh);
    return study;
}
}
