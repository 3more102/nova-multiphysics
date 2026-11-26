#include "test_framework.hpp"
#include "nova/physics/maxwell_cavity_3d.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

using namespace nova::physics;
using namespace nova::linalg;

NOVA_TEST(generalized_eigen_known_spectrum_mass_orthogonality) {
    const auto k=CsrMatrix::from_triplets(2,2,{{0,0,4},{0,1,2},{1,0,2},{1,1,4}});
    const auto m=CsrMatrix::from_triplets(2,2,{{0,0,2},{1,1,2}});
    const auto r=solve_symmetric_generalized_eigen(k,m);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE_NEAR(r.eigenpairs[0].eigenvalue,1.0,1e-12);
    NOVA_REQUIRE_NEAR(r.eigenpairs[1].eigenvalue,3.0,1e-12);
    for (std::size_t i=0;i<2;++i) {
        NOVA_REQUIRE(r.eigenpairs[i].relative_residual<1e-12);
        for (std::size_t j=0;j<2;++j) {
            const auto y=m.multiply(r.eigenpairs[j].eigenvector);
            double dot=0;
            for (std::size_t d=0;d<2;++d) dot+=y[d]*r.eigenpairs[i].eigenvector[d];
            NOVA_REQUIRE_NEAR(dot,i==j?1.0:0.0,1e-12);
        }
    }
}
NOVA_TEST(generalized_eigen_rejects_non_spd_and_asymmetry) {
    const auto eye=CsrMatrix::from_triplets(2,2,{{0,0,1},{1,1,1}});
    const auto indefinite=CsrMatrix::from_triplets(2,2,{{0,0,1},{1,1,-1}});
    const auto asym=CsrMatrix::from_triplets(2,2,{{0,1,1}});
    bool rejected=false;
    try {(void)solve_symmetric_generalized_eigen(eye,indefinite);}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    try {(void)solve_symmetric_generalized_eigen(asym,eye);}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
}
NOVA_TEST(maxwell_cavity_pec_cube_nullspace_and_eigenfrequencies) {
    constexpr double pi=std::numbers::pi;
    double previous=1.0;
    for (std::size_t n:{2U,3U,4U}) {
        const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,n,n,n).mesh;
        const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
        const auto a=assemble_maxwell_cavity_3d(mesh,materials);
        const auto r=solve_maxwell_cavity_3d(a,3);
        NOVA_REQUIRE(r.converged);
        NOVA_REQUIRE(r.nullspace_dimension==(n-1)*(n-1)*(n-1));
        NOVA_REQUIRE(r.modes.size()==3);
        const double reference=n==2?17.06363423:(n==3?18.43103487:18.96183604);
        NOVA_REQUIRE_NEAR(r.modes[0].vacuum_wavenumber_squared,reference,2e-7);
        const double expected_second=n==2?19.64300762:(n==3?20.02505572:19.94375703);
        for (std::size_t i=1;i<3;++i)
            NOVA_REQUIRE_NEAR(r.modes[i].vacuum_wavenumber_squared,expected_second,2e-7);
        const double error=std::abs(std::sqrt(r.modes[0].vacuum_wavenumber_squared/(2*pi*pi))-1);
        NOVA_REQUIRE(error<previous);
        previous=error;
        for (const auto& mode:r.modes) {
            NOVA_REQUIRE(mode.relative_residual<1e-10);
            for (std::size_t e=0;e<a.edges.size();++e)
                if (!std::binary_search(a.free_edges.begin(),a.free_edges.end(),e))
                    NOVA_REQUIRE_NEAR(mode.edge_field[e],0.0,0.0);
        }
        std::cout<<"PEC cube: cells="<<mesh.elements.size()<<" free="<<a.free_edges.size()
                 <<" nullity="<<r.nullspace_dimension<<" lambda1="
                 <<r.modes[0].vacuum_wavenumber_squared<<" frequency error="<<error<<'\n';
    }
}
NOVA_TEST(maxwell_cavity_dielectric_frequency_scaling_and_loss_rejection) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
    const auto a=solve_maxwell_cavity_3d(assemble_maxwell_cavity_3d(mesh,materials),1);
    for (auto& m:materials) m.relative_permittivity=4.0;
    const auto b=solve_maxwell_cavity_3d(assemble_maxwell_cavity_3d(mesh,materials),1);
    NOVA_REQUIRE(a.converged && b.converged);
    NOVA_REQUIRE_NEAR(b.modes[0].frequency_hz/a.modes[0].frequency_hz,0.5,1e-12);
    materials[0].relative_permittivity={4.0,0.1};
    bool rejected=false;
    try {(void)assemble_maxwell_cavity_3d(mesh,materials);}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
}


NOVA_TEST(maxwell_cavity_participation_energy_and_local_material_derivatives) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
    const auto a=assemble_maxwell_cavity_3d(mesh,materials);
    const auto modes=solve_maxwell_cavity_3d(a,2);
    NOVA_REQUIRE(modes.converged && modes.modes.size()==2);
    const auto& mode=modes.modes[0];
    NOVA_REQUIRE(modes.modes[1].vacuum_wavenumber_squared >
                 mode.vacuum_wavenumber_squared*1.01);
    const auto p=maxwell_cavity_participation_3d(mesh,materials,a,mode.edge_field);
    NOVA_REQUIRE(p.electric_fraction.size()==mesh.elements.size());
    NOVA_REQUIRE(p.magnetic_fraction.size()==mesh.elements.size());
    double se=0,sm=0;
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        NOVA_REQUIRE(p.electric_fraction[c]>=0.0 && p.magnetic_fraction[c]>=0.0);
        se+=p.electric_fraction[c]; sm+=p.magnetic_fraction[c];
    }
    NOVA_REQUIRE_NEAR(se,1.0,1e-12);
    NOVA_REQUIRE_NEAR(sm,1.0,1e-12);
    NOVA_REQUIRE_NEAR(p.mass_norm_squared,1.0,1e-11);
    NOVA_REQUIRE_NEAR(p.curl_energy/p.mass_norm_squared,
                      mode.vacuum_wavenumber_squared,1e-9);
    auto scaled=mode.edge_field;
    for (auto& v:scaled) v*=3;
    const auto q=maxwell_cavity_participation_3d(mesh,materials,a,scaled);
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        NOVA_REQUIRE_NEAR(q.electric_fraction[c],p.electric_fraction[c],1e-12);
        NOVA_REQUIRE_NEAR(q.magnetic_fraction[c],p.magnetic_fraction[c],1e-12);
    }
    constexpr double h=1e-4;
    for (bool electric:{true,false}) {
        const auto& fractions=electric?p.electric_fraction:p.magnetic_fraction;
        const auto c=static_cast<std::size_t>(
            std::max_element(fractions.begin(),fractions.end())-fractions.begin());
        auto frequency=[&](double delta) {
            auto modified=materials;
            if (electric) modified[c].relative_permittivity=std::exp(delta);
            else modified[c].relative_permeability=std::exp(delta);
            const auto result=solve_maxwell_cavity_3d(
                assemble_maxwell_cavity_3d(mesh,modified),1);
            NOVA_REQUIRE(result.converged && result.modes.size()==1);
            return result.modes[0].frequency_hz;
        };
        const double derivative=std::log(frequency(h)/frequency(-h))/(2*h);
        NOVA_REQUIRE_NEAR(derivative,-0.5*fractions[c],1e-7);
    }
    bool rejected=false;
    try {
        (void)maxwell_cavity_participation_3d(
            mesh,materials,a,std::vector<double>(a.edges.size()));
    } catch(const std::invalid_argument&) {rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    auto invalid=mode.edge_field;
    for (std::size_t e=0;e<a.edges.size();++e)
        if (!std::binary_search(a.free_edges.begin(),a.free_edges.end(),e)) {
            invalid[e]=1.0; break;
        }
    try {(void)maxwell_cavity_participation_3d(mesh,materials,a,invalid);}
    catch(const std::invalid_argument&) {rejected=true;}
    NOVA_REQUIRE(rejected);
}
