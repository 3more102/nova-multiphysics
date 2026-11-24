#include "nova/physics/package_3d.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc,char** argv) {
    using namespace nova;
    using Clock=std::chrono::steady_clock;
    try {
        if(argc!=3) throw std::invalid_argument("Usage: benchmark_package_3d subdivisions assemble|solve");
        std::size_t subdivisions{};
        const std::string_view arg=argv[1], mode=argv[2];
        const auto parsed=std::from_chars(arg.data(),arg.data()+arg.size(),subdivisions);
        if (parsed.ec!=std::errc{} || parsed.ptr!=arg.data()+arg.size() || subdivisions<1 ||
            (mode!="assemble" && mode!="solve")) throw std::invalid_argument("Invalid benchmark parameters.");
        const auto start=Clock::now();
        physics::PackageModel3D model;
        model.mesh=fem::make_box_mesh_3d(0,1,0,1,0,1,subdivisions,subdivisions,subdivisions).mesh;
        model.materials.push_back({300,{2,-.04},{4,.007},{2,.01},{100,0},.25,.001});
        model.cell_material.resize(model.mesh.elements.size());
        physics::PackageState3D initial;
        const auto n=model.mesh.nodes.size();
        initial.temperature.assign(n,300); initial.voltage.resize(n);
        for (std::size_t i=0;i<n;++i) {
            const auto p=model.mesh.nodes[i];
            initial.voltage[i]=2*p.x; // exact initial algebraic equilibrium
            if (p.x==0 || p.x==1) model.voltage.push_back({i,[p](double){return 2*p.x;}});
            if (p.x==0) model.supports.push_back({i,0,0});
            if (p.y==0) model.supports.push_back({i,1,0});
            if (p.z==0) model.supports.push_back({i,2,0});
        }
        physics::PreparedPackage3D study(std::move(model));
        const auto prepared=Clock::now();
        double residual=0, max_error=0;
        std::size_t nnz=0, newton=0, cg=0;
        {
            auto point=initial.temperature;
            point.insert(point.end(),initial.voltage.begin(),initial.voltage.end());
            auto system=study.assemble(point,initial,.01,time::implicit_formula(time::ImplicitScheme::BackwardEuler,.01),true);
            nnz=system.jacobian.non_zeros();
            if (std::abs(system.energy_defect_watts+8)>1e-7) throw std::runtime_error("Assembled source is incorrect.");
        }
        const auto assembled=Clock::now();
        if (mode=="solve") {
            physics::PackageOptions3D o;
            o.scheme=time::ImplicitScheme::BackwardEuler;
            o.initial_dt=o.minimum_dt=o.maximum_dt=.01;
            o.newton.gmres.solver.max_iterations=1500;
            o.newton.residual_tolerance=5e-9;
            study.advance(initial,.01,o,[&](const auto& s,const auto& out){
                residual=out.energy_defect_watts; newton=out.iterations; cg=out.mechanics.linear_solver.iterations;
                const double exact=300+std::log1p(.002)/.05;
                for (double T:s.temperature) max_error=std::max(max_error,std::abs(T-exact));
                if (max_error>5e-5 || std::abs(residual)>2e-6 ||
                    std::abs(out.transport.terminal_power-out.transport.joule_power)>2e-6)
                    throw std::runtime_error("Coupled benchmark failed conservation/analytic check.");
            });
        }
        const auto elapsed=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
        std::cout<<std::setprecision(12)<<"nodes="<<n<<" total_field_dofs="<<5*n
                 <<" transport_rows="<<2*n<<" tets="<<study.model().mesh.elements.size()<<" nnz="<<nnz
                 <<" setup_s="<<elapsed(start,prepared)<<" assembly_s="<<elapsed(prepared,assembled)
                 <<" advance_s="<<elapsed(assembled,Clock::now())<<" newton="<<newton<<" cg="<<cg
                 <<" defect_W="<<residual<<" temperature_error_K="<<max_error<<" mode="<<mode<<'\n';
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
