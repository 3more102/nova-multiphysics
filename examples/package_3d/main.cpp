#include "nova/io/package_study_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_series.hpp"
#include "nova/fem/tetrahedron_h1.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace fs=std::filesystem;
namespace {
std::ofstream output(const fs::path& path) {
    std::ofstream out(path,std::ios::binary);
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<std::setprecision(17);
    return out;
}
std::string name(std::size_t step) {
    std::ostringstream out; out<<"step_"<<std::setfill('0')<<std::setw(6)<<step; return out.str();
}
}
int main(int argc,char** argv) {
    using namespace nova;
    try {
        if (argc<3) throw std::invalid_argument("Usage: nova_package_3d model.nova new_output_dir [--restart state.chk] [--until seconds]");
        auto study=io::read_package_study_3d(argv[1]);
        fs::path directory=argv[2], restart;
        double end=study.end_time;
        for (int i=3;i<argc;i+=2) {
            if (i+1>=argc) throw std::invalid_argument("Missing command-line value.");
            const std::string key=argv[i], value=argv[i+1];
            if (key=="--restart" && restart.empty()) restart=value;
            else if (key=="--until") {
                const auto parsed=std::from_chars(value.data(),value.data()+value.size(),end);
                if (parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size() || !std::isfinite(end) || end>study.end_time)
                    throw std::invalid_argument("Invalid --until time.");
            } else throw std::invalid_argument("Unexpected option "+key);
        }
        physics::PreparedPackage3D problem(std::move(study.model));
        const auto& mesh=problem.model().mesh;
        physics::PackageState3D state;
        if (restart.empty()) state=problem.initialize(std::vector<double>(mesh.nodes.size(),study.initial_kelvin));
        else {
            std::ifstream in(restart,std::ios::binary);
            state=io::read_package_checkpoint_3d(in,mesh,study.model_id);
        }
        if (!(end>state.time)) throw std::invalid_argument("End time must exceed the input state's time.");
        // Claim a run directory atomically. Independent runs never share or replace evidence.
        if (!fs::create_directories(directory)) throw std::invalid_argument("Output directory already exists.");
        auto history=output(directory/"history.csv");
        auto currents=output(directory/"terminals.csv");
        history<<"time_s,dt_s,rejected,newton,scaled_residual,error_ratio,Tmin_K,Tmax_K,Tmean_K,enthalpy_J,joule_W,source_W,cooling_W,dirichlet_W,heat_defect_W,terminal_W,net_current_A,max_u_m,max_von_mises_Pa";
        for (const auto& p:study.probes) history<<','<<p.name<<"_K,"<<p.name<<"_V";
        history<<'\n'; currents<<"time_s,terminal_tag,current_into_domain_A\n";
        std::vector<double> volume(mesh.elements.size());
        for (std::size_t e=0;e<volume.size();++e)
            volume[e]=std::abs(fem::evaluate_tetrahedron_p1_physical_gradients_3d(mesh,e).signed_jacobian_determinant)/6;
        const double total_volume=std::accumulate(volume.begin(),volume.end(),0.);
        std::vector<io::VtkTimeSeriesEntry> frames;
        std::size_t count=0,rejected=0;
        double max_defect=0,max_power_error=0,last_max=0;
        const auto start=std::chrono::steady_clock::now();
        problem.advance(state,end,study.options,[&](const auto& s,const auto& step) {
            ++count; rejected+=step.rejected_before;
            const auto& f=step.transport;
            const auto [low,high]=std::minmax_element(s.temperature.begin(),s.temperature.end());
            last_max=*high;
            max_defect=std::max(max_defect,std::abs(step.energy_defect_watts));
            max_power_error=std::max(max_power_error,std::abs(f.terminal_power-f.joule_power));
            double mean=0,max_u=0,max_vm=0;
            for (std::size_t e=0;e<volume.size();++e) {
                for (auto i:mesh.elements[e]) mean+=.25*volume[e]*s.temperature[i];
                max_vm=std::max(max_vm,physics::von_mises_stress_3d(step.mechanics.cell_stress[e]));
            }
            for (const auto& u:step.mechanics.displacement) max_u=std::max(max_u,std::hypot(u[0],u[1],u[2]));
            history<<s.time<<','<<s.time-s.previous_time<<','<<step.rejected_before<<','<<step.iterations<<','
                   <<step.residual<<','<<step.normalized_error<<','<<*low<<','<<*high<<','<<mean/total_volume<<','
                   <<f.enthalpy<<','<<f.joule_power<<','<<f.source_power<<','<<f.cooling_power<<','
                   <<step.dirichlet_heat_watts<<','<<step.energy_defect_watts<<','<<f.terminal_power<<','
                   <<f.net_terminal_current<<','<<max_u<<','<<max_vm;
            for (const auto& p:study.probes) {
                double T=0,V=0;
                for (unsigned j=0;j<4;++j) {
                    const auto node=mesh.elements[p.cell][j];
                    T+=p.weights[j]*s.temperature[node]; V+=p.weights[j]*s.voltage[node];
                }
                history<<','<<T<<',';
                if (problem.model().materials[problem.model().cell_material[p.cell]].electrical.value>0) history<<V;
            }
            history<<'\n'; history.flush();
            std::map<int,double> terminal;
            for (std::size_t j=0;j<f.terminal_current.size();++j) terminal[study.terminal_tags[j]]+=f.terminal_current[j];
            for (auto [tag,I]:terminal) currents<<s.time<<','<<tag<<','<<I<<'\n';
            currents.flush();
            const std::string prefix=name(count);
            if (study.frame_stride && (count%study.frame_stride==0 || s.time==end)) {
                const auto file=prefix+".vtk";
                auto vtk=output(directory/file); io::write_legacy_vtk_package_3d(vtk,mesh,study.material_tags,s,step);
                vtk.close(); frames.push_back({s.time,file});
            }
            if (s.time==end || std::binary_search(study.options.events.begin(),study.options.events.end(),s.time)) {
                const auto staging=directory/(prefix+".tmp"), file=directory/(prefix+".chk");
                auto chk=output(staging); io::write_package_checkpoint_3d(chk,mesh,study.model_id,s); chk.close();
                std::ifstream in(staging,std::ios::binary);
                const auto restored=io::read_package_checkpoint_3d(in,mesh,study.model_id);
                if (restored.temperature!=s.temperature || restored.voltage!=s.voltage ||
                    restored.previous_temperature!=s.previous_temperature || restored.previous_time!=s.previous_time ||
                    restored.time!=s.time || restored.next_dt!=s.next_dt) throw std::runtime_error("Checkpoint read-back mismatch.");
                in.close(); fs::rename(staging,file);
            }
        });
        if (!frames.empty()) { auto out=output(directory/"package.pvd"); io::write_pvd_collection(out,frames); out.close(); }
        auto summary=output(directory/"summary.json");
        summary<<"{\n  \"model_id\": \""<<study.model_id<<"\",\n  \"nodes\": "<<mesh.nodes.size()
               <<",\n  \"tetrahedra\": "<<mesh.elements.size()<<",\n  \"accepted\": "<<count
               <<",\n  \"rejected\": "<<rejected<<",\n  \"final_time_s\": "<<state.time
               <<",\n  \"final_max_temperature_K\": "<<last_max<<",\n  \"max_heat_defect_W\": "<<max_defect
               <<",\n  \"max_electrical_power_defect_W\": "<<max_power_error<<",\n  \"advance_and_export_s\": "
               <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";
        summary.close(); history.close(); currents.close();
        std::cout<<std::setprecision(12)<<"Completed t="<<state.time<<" s, accepted="<<count<<", rejected="<<rejected
                 <<", Tmax="<<last_max<<" K, max heat defect="<<max_defect<<" W\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
