#include "package_transport.hpp"
#include "nova/fem/simplex_quadrature.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {
using Vec = ElasticVector3D;
double dot(const Vec& a, const Vec& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec transform(const fem::SymmetricTensor3D& D, const Vec& g) {
    return {D.xx*g[0]+D.xy*g[1]+D.xz*g[2], D.xy*g[0]+D.yy*g[1]+D.yz*g[2],
            D.xz*g[0]+D.yz*g[1]+D.zz*g[2]};
}
double property(ExponentialProperty3D p, double delta) {
    if (p.value == 0) return 0;
    const double value = p.value * std::exp(p.slope*delta);
    if (!(value > 0) || !std::isfinite(value)) throw std::domain_error("Package material outside its finite positive range.");
    return value;
}
// Evaluate H(right)-H(left) accurately even for small temperature increments.
double enthalpy_increment(ExponentialProperty3D c, double ref, double right, double left) {
    if (c.slope == 0) return c.value*(right-left);
    return property(c,left-ref)*(std::expm1(c.slope*(right-left))/c.slope);
}
double signal(const TimeSignal3D& f, double t) {
    const double v = f(t);
    if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite package boundary or source signal.");
    return v;
}
void require_field(std::span<const double> v, std::size_t n, bool kelvin) {
    if (v.size()!=n) throw std::invalid_argument("Package field size disagrees with mesh.");
    for (double x : v) if (!std::isfinite(x) || (kelvin && x <= 0))
        throw std::domain_error("Expected finite field values and positive absolute temperatures.");
}
double area(const fem::Point3D& a, const fem::Point3D& b, const fem::Point3D& c) {
    const Vec u{b.x-a.x,b.y-a.y,b.z-a.z}, v{c.x-a.x,c.y-a.y,c.z-a.z};
    return 0.5*std::hypot(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]);
}
auto sorted(fem::BoundaryFace3D f) { std::sort(f.begin(),f.end()); return f; }
template<class Sample> void face_rule(bool radiation, const Sample& sample) {
    if (!radiation) {
        for (const auto& q : fem::triangle_face_degree_two_quadrature) sample(q.shape,2*q.weight);
        return;
    }
    // Tensor Gauss rule under a Duffy map. Integrates all face polynomials through degree five.
    constexpr double u[]{.06943184420297371,.33000947820757187,.6699905217924281,.9305681557970262};
    constexpr double wu[]{.17392742256872693,.32607257743127307,.32607257743127307,.17392742256872693};
    constexpr double v[]{.11270166537925831,.5,.8872983346207417};
    constexpr double wv[]{5./18,4./9,5./18};
    for (unsigned i=0;i<4;++i) for (unsigned j=0;j<3;++j)
        sample(std::array<double,3>{1-u[i],u[i]*(1-v[j]),u[i]*v[j]},2*wu[i]*wv[j]*u[i]);
}
}

PreparedPackage3D::Impl::Impl(PackageModel3D input) : model(std::move(input)) {
    const auto& mesh=model.mesh;
    const auto topology=fem::analyze_tetrahedral_topology_3d(mesh);
    const auto n=mesh.nodes.size();
    if (mesh.elements.empty() || model.cell_material.size()!=mesh.elements.size() ||
        (!model.heat_sources.empty() && model.heat_sources.size()!=model.materials.size()))
        throw std::invalid_argument("Package mesh, cell materials or heat-source table is incomplete.");
    for (const auto& m : model.materials) {
        if (!std::isfinite(m.reference_temperature) || m.reference_temperature<=0 ||
            !std::isfinite(m.expansion) || !(m.poisson > -1 && m.poisson < 0.5) ||
            m.electrical.value<0 || !(m.thermal.value>0 && m.young.value>0 && m.volumetric_capacity.value>0))
            throw std::invalid_argument("Invalid reference package material.");
        if (m.electrical_shape) fem::validate_conductivity_tensor_3d(*m.electrical_shape);
        if (m.thermal_shape) fem::validate_conductivity_tensor_3d(*m.thermal_shape);
        for (const auto p : {m.electrical,m.thermal,m.volumetric_capacity,m.young})
            if (!std::isfinite(p.value) || !std::isfinite(p.slope))
                throw std::invalid_argument("Nonfinite package material parameter.");
    }
    active.assign(n,0); fixed.assign(2*n,0);
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(),parent.end(),0U);
    const auto root=[&](std::size_t i) {
        while (parent[i]!=i) { parent[i]=parent[parent[i]]; i=parent[i]; }
        return i;
    };
    cells.reserve(mesh.elements.size());
    for (std::size_t e=0;e<mesh.elements.size();++e) {
        const auto& m=model.materials.at(model.cell_material[e]);
        const auto basis=fem::evaluate_tetrahedron_p1_physical_gradients_3d(mesh,e);
        cells.push_back({basis.gradients,std::abs(basis.signed_jacobian_determinant)/6});
        expansion.push_back({m.expansion,m.reference_temperature});
        variable_modulus |= m.young.slope!=0;
        if (m.electrical.value>0) {
            const auto& nodes=mesh.elements[e];
            for (auto i:nodes) { active[i]=1; parent[root(i)]=root(nodes[0]); }
        }
    }
    const auto bind=[&](const auto& list,std::size_t offset) {
        for (const auto& bc:list) {
            if (bc.node>=n || !bc.value || fixed[offset+bc.node] || (offset && !active[bc.node]))
                throw std::invalid_argument("Invalid, inactive or duplicate package Dirichlet node.");
            fixed[offset+bc.node]=1;
        }
    };
    bind(model.temperature,0); bind(model.voltage,n);
    std::vector<bool> grounded(n);
    for (const auto& bc:model.voltage) grounded[root(bc.node)]=true;
    for (std::size_t i=0;i<n;++i) {
        if (active[i] && !grounded[root(i)]) throw std::invalid_argument("Each conducting component requires a voltage reference.");
        if (!active[i]) fixed[n+i]=1;
    }
    std::vector<std::vector<std::size_t>> neighbors(n);
    for (const auto& tet:mesh.elements) for (auto i:tet)
        neighbors[i].insert(neighbors[i].end(),tet.begin(),tet.end());
    for (auto& nodes:neighbors) {
        std::sort(nodes.begin(),nodes.end());
        nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
    }
    std::vector<std::size_t> offsets{0}, columns;
    offsets.reserve(2*n+1);
    for (std::size_t row=0;row<2*n;++row) {
        if (fixed[row]) columns.push_back(row);
        else {
            const auto& nodes=neighbors[row%n];
            columns.insert(columns.end(),nodes.begin(),nodes.end());
            for (auto i:nodes) columns.push_back(n+i);
        }
        offsets.push_back(columns.size());
    }
    pattern=linalg::CsrAssemblyPattern(2*n,std::move(offsets),std::move(columns));
    std::set<fem::BoundaryFace3D> exterior;
    for (auto f:topology.exterior_faces) exterior.insert(sorted(f));
    for (const auto& bc:model.cooling) {
        if (!exterior.erase(sorted(bc.face)) || !std::isfinite(bc.convection) || bc.convection<0 ||
            !(bc.emissivity>=0 && bc.emissivity<=1) || !bc.ambient_kelvin)
            throw std::invalid_argument("Invalid, repeated or interior package cooling face.");
        areas.push_back(area(mesh.nodes[bc.face[0]],mesh.nodes[bc.face[1]],mesh.nodes[bc.face[2]]));
    }
}

void PreparedPackage3D::Impl::integrate(std::span<const double> temp, std::span<const double> volt,
    double t, const PackageState3D* old, time::ImplicitFormula f, bool tangent,
    bool electric_rows, Evaluation& out, TransportFields3D* fields) const {
    const auto n=model.mesh.nodes.size();
    const auto add=[&](std::size_t row,std::size_t col,double value) {
        if (!fixed[row] && value!=0) out.entries[pattern.slot(row,col)]+=value;
    };
    std::vector<double> sources(model.materials.size());
    for (std::size_t i=0;i<model.heat_sources.size();++i)
        if (model.heat_sources[i]) sources[i]=signal(model.heat_sources[i],t);
    for (std::size_t e=0;e<cells.size();++e) {
        const auto& nodes=model.mesh.elements[e];
        const auto& cell=cells[e];
        const auto& grad=cell.gradient;
        const auto& mat=model.materials[model.cell_material[e]];
        Vec gt{},gv{};
        double mean=0;
        for (unsigned b=0;b<4;++b) {
            mean+=0.25*temp[nodes[b]];
            for (unsigned d=0;d<3;++d) {
                gt[d]+=(temp[nodes[b]]-temp[nodes[0]])*grad[b][d];
                gv[d]+=(volt[nodes[b]]-volt[nodes[0]])*grad[b][d];
            }
        }
        const double delta=mean-mat.reference_temperature;
        const double sigma=property(mat.electrical,delta), k=property(mat.thermal,delta);
        const Vec flux_v=mat.electrical_shape ? transform(*mat.electrical_shape,gv) : gv;
        const Vec flux_t=mat.thermal_shape ? transform(*mat.thermal_shape,gt) : gt;
        const double joule=sigma*dot(gv,flux_v), Q=sources[model.cell_material[e]];
        const double w=cell.volume;
        std::array<double,4> storage{}, capacity{}, dgV{}, dgT{};
        std::array<Vec,4> thermal_basis{}, electrical_basis{};
        for (unsigned b=0;b<4;++b) {
            const auto j=nodes[b];
            dgV[b]=dot(grad[b],flux_v); dgT[b]=dot(grad[b],flux_t);
            if (tangent && mat.thermal_shape) thermal_basis[b]=transform(*mat.thermal_shape,grad[b]);
            if (tangent && mat.electrical_shape) electrical_basis[b]=transform(*mat.electrical_shape,grad[b]);
            if (old) {
                storage[b]=f.a*enthalpy_increment(mat.volumetric_capacity,mat.reference_temperature,temp[j],old->temperature[j]);
                if (f.b) storage[b]-=f.b*enthalpy_increment(mat.volumetric_capacity,mat.reference_temperature,
                                                           old->temperature[j],old->previous_temperature[j]);
                if (tangent) capacity[b]=f.a*property(mat.volumetric_capacity,temp[j]-mat.reference_temperature);
            }
            if (fields) fields->enthalpy+=0.25*w*enthalpy_increment(mat.volumetric_capacity,mat.reference_temperature,
                                                                   temp[j],mat.reference_temperature);
        }
        const double sum_storage=std::accumulate(storage.begin(),storage.end(),0.);
        for (unsigned a=0;a<4;++a) {
            const auto i=nodes[a];
            out.residual[i]+=w*((sum_storage+storage[a])/20 + f.theta*(k*dgT[a]-(joule+Q)/4));
            if (electric_rows) out.residual[n+i]+=w*sigma*dgV[a];
            if (!tangent) continue;
            for (unsigned b=0;b<4;++b) {
                const auto j=nodes[b];
                const double G=dot(grad[a],grad[b]);
                const double Gt=mat.thermal_shape ? dot(grad[a],thermal_basis[b]) : G;
                const double Gv=mat.electrical_shape ? dot(grad[a],electrical_basis[b]) : G;
                add(i,j,w*((a==b ? .1 : .05)*capacity[b]+f.theta*
                    (k*Gt+k*mat.thermal.slope*dgT[a]/4-joule*mat.electrical.slope/16)));
                add(i,n+j,-w*f.theta*sigma*dgV[b]/2);
                add(n+i,j,w*sigma*mat.electrical.slope*dgV[a]/4);
                add(n+i,n+j,w*sigma*Gv);
            }
        }
        if (fields) {
            fields->heat_flux.push_back({{gt[0],gt[1],gt[2]},{-k*flux_t[0],-k*flux_t[1],-k*flux_t[2]}});
            fields->current_density.push_back({-sigma*flux_v[0],-sigma*flux_v[1],-sigma*flux_v[2]});
            fields->joule_density.push_back(joule);
            fields->conductivity.push_back(sigma);
            fields->young_modulus.push_back(property(mat.young,delta));
            fields->joule_power+=w*joule;
            fields->source_power+=w*Q;
        }
    }
    constexpr double stefan_boltzmann=5.670374419e-8;
    for (std::size_t i=0;i<model.cooling.size();++i) {
        const auto& bc=model.cooling[i];
        const double ambient=signal(bc.ambient_kelvin,t);
        if (!(ambient>0)) throw std::invalid_argument("Radiation and convection require absolute ambient temperature.");
        const double external=bc.outward_flux ? signal(bc.outward_flux,t) : 0;
        face_rule(bc.emissivity>0,[&](const auto& N,double weight) {
            double T=0;
            for (unsigned a=0;a<3;++a) T+=N[a]*temp[bc.face[a]];
            const double radiation=bc.emissivity*stefan_boltzmann;
            const double flux=external+(T-ambient)*(bc.convection+radiation*(T+ambient)*(T*T+ambient*ambient));
            const double derivative=bc.convection+4*radiation*T*T*T;
            const double measure=areas[i]*weight;
            if (fields) fields->cooling_power+=measure*flux;
            for (unsigned a=0;a<3;++a) {
                const auto row=bc.face[a];
                out.residual[row]+=f.theta*measure*N[a]*flux;
                if (tangent) for (unsigned b=0;b<3;++b)
                    add(row,bc.face[b],f.theta*measure*N[a]*N[b]*derivative);
            }
        });
    }
}

std::vector<double> PreparedPackage3D::Impl::old_operator(const PackageState3D& old,double theta) const {
    if (theta==1) return {};
    Evaluation past;
    past.residual.resize(2*model.mesh.nodes.size());
    integrate(old.temperature,old.voltage,old.time,nullptr,{0,0,1-theta},false,false,past);
    past.residual.resize(model.mesh.nodes.size());
    return past.residual;
}
TransportAssembly3D PreparedPackage3D::Impl::system(std::span<const double> x,const PackageState3D& old,
    double t,time::ImplicitFormula f,std::span<const double> past,bool tangent) const {
    const auto n=model.mesh.nodes.size();
    Evaluation out;
    out.residual.resize(2*n);
    if (tangent) out.entries.resize(pattern.slots());
    std::copy(past.begin(),past.end(),out.residual.begin());
    integrate(x.first(n),x.subspan(n),t,&old,f,tangent,true,out);
    TransportAssembly3D result;
    long double defect=0, reaction=0;
    for (std::size_t i=0;i<n;++i) {
        if (fixed[i]) reaction+=out.residual[i]; else defect+=out.residual[i];
        if (!active[i]) out.residual[n+i]=x[n+i];
    }
    result.energy_defect_watts=static_cast<double>(defect);
    result.dirichlet_heat_watts=static_cast<double>(reaction);
    for (const auto& bc:model.temperature) out.residual[bc.node]=x[bc.node]-signal(bc.value,t);
    for (const auto& bc:model.voltage) out.residual[n+bc.node]=x[n+bc.node]-signal(bc.value,t);
    if (tangent) {
        for (std::size_t i=0;i<2*n;++i) if (fixed[i]) out.entries[pattern.slot(i,i)]=1;
        result.jacobian=pattern.matrix(out.entries);
    }
    result.residual=std::move(out.residual);
    return result;
}
PreparedPackage3D::PreparedPackage3D(PackageModel3D m) : impl_(std::make_unique<Impl>(std::move(m))) {}
PreparedPackage3D::~PreparedPackage3D()=default;
PreparedPackage3D::PreparedPackage3D(PreparedPackage3D&&) noexcept=default;
PreparedPackage3D& PreparedPackage3D::operator=(PreparedPackage3D&&) noexcept=default;
const PackageModel3D& PreparedPackage3D::model() const noexcept { return impl_->model; }

TransportAssembly3D PreparedPackage3D::assemble(std::span<const double> x,const PackageState3D& old,
    double t,time::ImplicitFormula f,bool tangent) const {
    const auto n=model().mesh.nodes.size();
    if (x.size()!=2*n || !std::isfinite(t) || !std::isfinite(old.time) || !(t>old.time) ||
        !(f.a>0) || !std::isfinite(f.a) || !std::isfinite(f.b) || f.b<0 || !(f.theta>0 && f.theta<=1))
        throw std::invalid_argument("Invalid package assembly dimensions or time formula.");
    require_field(x.first(n),n,true); require_field(x.subspan(n),n,false);
    require_field(old.temperature,n,true); require_field(old.voltage,n,false);
    if (f.b) require_field(old.previous_temperature,n,true);
    return impl_->system(x,old,t,f,impl_->old_operator(old,f.theta),tangent);
}
TransportFields3D PreparedPackage3D::recover(std::span<const double> T,std::span<const double> V,double t) const {
    const auto n=model().mesh.nodes.size();
    require_field(T,n,true); require_field(V,n,false);
    if (!std::isfinite(t)) throw std::invalid_argument("Nonfinite recovery time.");
    TransportFields3D fields;
    fields.electrically_active=impl_->active;
    Impl::Evaluation out;
    out.residual.resize(2*n);
    impl_->integrate(T,V,t,nullptr,{0,0,1},false,true,out,&fields);
    for (const auto& bc:model().voltage) {
        const double current=out.residual[n+bc.node];
        fields.terminal_current.push_back(current);
        fields.net_terminal_current+=current;
        fields.terminal_power+=V[bc.node]*current;
    }
    return fields;
}
LinearElasticity3DResult PreparedPackage3D::Impl::equilibrate(std::span<const double> T,
    const linalg::SolverOptions& options) {
    if (!elasticity || variable_modulus) {
        std::vector<IsotropicElasticMaterial3D> laws;
        laws.reserve(cells.size());
        for (std::size_t e=0;e<cells.size();++e) {
            const auto& m=model.materials[model.cell_material[e]];
            double delta=0;
            for (auto i:model.mesh.elements[e]) delta+=.25*(T[i]-m.reference_temperature);
            laws.push_back({property(m.young,delta),m.poisson});
        }
        elasticity=std::make_unique<PreparedElasticity3D>(model.mesh,laws,model.clamps,
            model.body_forces,model.tractions,model.supports);
    }
    return elasticity->solve_temperature(T,expansion,options);
}
}
