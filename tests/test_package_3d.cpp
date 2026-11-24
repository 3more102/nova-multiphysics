#include "nova/physics/package_3d.hpp"
#include "test_framework.hpp"
#include <algorithm>
#include <cmath>

namespace {
using namespace nova;
using namespace nova::physics;
PackageModel3D cube(std::size_t size=2) {
    PackageModel3D m;
    m.mesh=fem::make_box_mesh_3d(0,1,0,1,0,1,size,1,1).mesh;
    m.materials.push_back({300,{2,-.04},{4,.007},{2,.01},{100,-.013},.25,.001});
    m.cell_material.resize(m.mesh.elements.size());
    for (std::size_t i=0;i<m.mesh.nodes.size();++i) {
        const auto p=m.mesh.nodes[i];
        if (p.x==0 || p.x==1) m.voltage.push_back({i,[p](double){ return 2*p.x; }});
        if (p.x==0) m.supports.push_back({i,0,0});
        if (p.y==0) m.supports.push_back({i,1,0});
        if (p.z==0) m.supports.push_back({i,2,0});
    }
    return m;
}
PackageOptions3D fixed(double dt,time::ImplicitScheme scheme) {
    PackageOptions3D o;
    o.initial_dt=dt; o.maximum_dt=dt; o.scheme=scheme;
    return o;
}
double exact(double t) { return 300+std::log1p(.2*t)/.05; }
template<class F> bool throws(F f) { try { f(); } catch(const std::exception&) { return true; } return false; }
}

NOVA_TEST(package_variable_step_bdf2_differentiates_quadratics_exactly) {
    for (double r:{.17,.5,1.,1.7,2.}) {
        const double h=.14*r, k=.14, t=.43;
        const auto f=time::implicit_formula(time::ImplicitScheme::BDF2,h,k);
        for (const auto& p:std::vector<std::pair<double,double>>{{1,0},{0,1},{1.7,-.4}}) {
            const auto y=[&](double s){ return p.first*s*s+p.second*s; };
            NOVA_REQUIRE_NEAR(f.a*(y(t)-y(t-h))-f.b*(y(t-h)-y(t-h-k)),2*p.first*t+p.second,2e-14);
        }
    }
}

NOVA_TEST(package_full_tangent_matches_independent_column_differences) {
    auto m=cube();
    m.temperature.push_back({0,[](double t){ return 300+t; }});
    for (auto face:fem::extract_boundary_faces_3d(m.mesh))
        if (std::all_of(face.begin(),face.end(),[&](auto i){ return m.mesh.nodes[i].y==1; }))
            m.cooling.push_back({face,8,.8,[](double){return 290.;},[](double t){return 2*t;}});
    PreparedPackage3D p(std::move(m));
    const auto n=p.model().mesh.nodes.size();
    auto old=p.initialize(std::vector<double>(n,300));
    old.previous_temperature.assign(n,299);
    std::vector<double> x(2*n);
    for (std::size_t i=0;i<n;++i) {
        x[i]=301+.2*static_cast<double>(i);
        x[n+i]=old.voltage[i]+.13*static_cast<double>(i);
    }
    for (const auto scheme:{time::ImplicitScheme::BDF2,time::ImplicitScheme::CrankNicolson}) {
        const auto f=time::implicit_formula(scheme,.2,.15);
        auto J=p.assemble(x,old,.2,f,true).jacobian;
        for (std::size_t j=0;j<2*n;++j) {
            const double h=j<n ? .0001 : .00001;
            auto plus=x,minus=x; plus[j]+=h; minus[j]-=h;
            const auto rp=p.assemble(plus,old,.2,f,false).residual;
            const auto rm=p.assemble(minus,old,.2,f,false).residual;
            for (std::size_t i=0;i<2*n;++i) {
                const double fd=(rp[i]-rm[i])/(2*h);
                NOVA_REQUIRE_NEAR(J.at(i,j),fd,1e-8*(1+std::abs(fd)));
            }
        }
    }
}
NOVA_TEST(package_nonlinear_joule_enthalpy_temporal_order_and_free_expansion) {
    PreparedPackage3D p(cube());
    for (const auto scheme:{time::ImplicitScheme::BackwardEuler,time::ImplicitScheme::CrankNicolson,time::ImplicitScheme::BDF2}) {
        double previous_error=0;
        for (double h:{.1,.05,.025}) {
            auto state=p.initialize(std::vector<double>(p.model().mesh.nodes.size(),300.));
            p.advance(state,.5,fixed(h,scheme),[&](const auto& s,const auto& out) {
                const auto& fields=out.transport;
                NOVA_REQUIRE_NEAR(fields.terminal_power,fields.joule_power,1e-11);
                NOVA_REQUIRE_NEAR(fields.net_terminal_current,0,1e-11);
                NOVA_REQUIRE_NEAR(out.energy_defect_watts,0,1e-8);
                NOVA_REQUIRE_NEAR(fields.conductivity[0],2*std::exp(-.04*(s.temperature[0]-300)),1e-11);
                NOVA_REQUIRE_NEAR(fields.young_modulus[0],100*std::exp(-.013*(s.temperature[0]-300)),1e-10);
                for (std::size_t i=0;i<s.voltage.size();++i) {
                    const auto point=p.model().mesh.nodes[i];
                    NOVA_REQUIRE_NEAR(s.voltage[i],2*point.x,1e-12);
                    NOVA_REQUIRE_NEAR(s.temperature[i],s.temperature[0],5e-10);
                    NOVA_REQUIRE_NEAR(out.mechanics.displacement[i][2],.001*(s.temperature[0]-300)*point.z,1e-11);
                }
                for (const auto& stress:out.mechanics.cell_stress) for (double value:stress)
                    NOVA_REQUIRE_NEAR(value,0,2e-9);
            });
            double error=std::abs(state.temperature[0]-exact(.5));
            if (previous_error) NOVA_REQUIRE(previous_error/error>(scheme==time::ImplicitScheme::BackwardEuler ? 1.9 : 3.4));
            previous_error=error;
        }
    }
}
NOVA_TEST(package_segregated_monolithic_equivalence_on_nonuniform_temperature) {
    auto m=cube(3);
    for (std::size_t i=0;i<m.mesh.nodes.size();++i) if (m.mesh.nodes[i].x==0)
        m.temperature.push_back({i,[](double t){ return 300+5*t; }});
    PreparedPackage3D p(std::move(m));
    auto a=p.initialize(std::vector<double>(p.model().mesh.nodes.size(),300.)), b=a;
    auto o=fixed(.1,time::ImplicitScheme::BDF2);
    p.advance(a,.3,o);
    o.strategy=TransportStrategy3D::Segregated; o.relaxation=.8;
    p.advance(b,.3,o);
    for (std::size_t i=0;i<a.voltage.size();++i) {
        NOVA_REQUIRE_NEAR(a.temperature[i],b.temperature[i],2e-9);
        NOVA_REQUIRE_NEAR(a.voltage[i],b.voltage[i],2e-10);
    }
}
NOVA_TEST(package_event_adaptive_rollback_and_accepted_restart) {
    PreparedPackage3D p(cube());
    auto initial=p.initialize(std::vector<double>(p.model().mesh.nodes.size(),300.));
    auto state=initial;
    auto o=fixed(.3,time::ImplicitScheme::BDF2);
    o.events={.17}; o.error_kelvin=.0004;
    bool event_seen=false; std::size_t rejections=0;
    PackageState3D at_event;
    p.advance(state,.5,o,[&](const auto& s,const auto& out){
        NOVA_REQUIRE(out.normalized_error<=1);
        rejections+=out.rejected_before;
        if (s.time==.17) { event_seen=true; at_event=s; }
    });
    NOVA_REQUIRE(event_seen && rejections>0);
    p.advance(at_event,.5,o);
    NOVA_REQUIRE(at_event.temperature==state.temperature && at_event.voltage==state.voltage);
    NOVA_REQUIRE_NEAR(state.temperature[0],exact(.5),.0005);
    o.minimum_dt=.3; o.maximum_attempts=1;
    NOVA_REQUIRE(throws([&]{ p.advance(initial,.5,o); }));
    NOVA_REQUIRE(initial.time==0 && initial.temperature[0]==300 && initial.previous_temperature.empty());
}
NOVA_TEST(package_exact_insulator_and_unanchored_conductor_validation) {
    auto m=cube();
    m.voltage.clear();
    NOVA_REQUIRE(throws([&]{PreparedPackage3D bad(m);}));
    m.materials[0].electrical.value=0;
    PreparedPackage3D p(std::move(m));
    auto s=p.initialize(std::vector<double>(p.model().mesh.nodes.size(),300));
    p.advance(s,.1,fixed(.1,time::ImplicitScheme::BackwardEuler),[](const auto&,const auto& out){
        NOVA_REQUIRE(out.transport.joule_power==0);
        for (auto active:out.transport.electrically_active) NOVA_REQUIRE(active==0);
    });
}
