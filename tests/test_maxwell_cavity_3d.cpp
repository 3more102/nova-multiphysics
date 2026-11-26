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

NOVA_TEST(maxwell_cavity_geometric_and_magnetic_similarity) {
    // For geometrically similar PEC cavities, K scales as 1/L and M as L.
    // With uniform relative permeability mu_r, K also scales as 1/mu_r.
    // Thus lambda scales as 1/(L^2 mu_r), and f as 1/(L sqrt(mu_r)).
    const auto unit_mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    const auto large_mesh=nova::fem::make_box_mesh_3d(0,2,0,2,0,2,2,2,2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(unit_mesh.elements.size());
    const auto baseline=solve_maxwell_cavity_3d(
        assemble_maxwell_cavity_3d(unit_mesh,materials),3);
    NOVA_REQUIRE(baseline.converged);
    NOVA_REQUIRE(baseline.modes.size()==3);
    const auto geometric=solve_maxwell_cavity_3d(
        assemble_maxwell_cavity_3d(large_mesh,materials),3);
    for (auto& material:materials) material.relative_permeability=4.0;
    const auto magnetic=solve_maxwell_cavity_3d(
        assemble_maxwell_cavity_3d(unit_mesh,materials),3);
    const auto combined=solve_maxwell_cavity_3d(
        assemble_maxwell_cavity_3d(large_mesh,materials),3);
    for (const auto* result:{&geometric,&magnetic,&combined}) {
        NOVA_REQUIRE(result->converged);
        NOVA_REQUIRE(result->modes.size()==baseline.modes.size());
        NOVA_REQUIRE(result->nullspace_dimension==baseline.nullspace_dimension);
        for (const auto& mode:result->modes)
            NOVA_REQUIRE(mode.relative_residual<1e-10);
    }
    for (std::size_t i=0;i<baseline.modes.size();++i) {
        const auto& ref=baseline.modes[i];
        NOVA_REQUIRE(ref.vacuum_wavenumber_squared>0.0);
        NOVA_REQUIRE(ref.frequency_hz>0.0);
        NOVA_REQUIRE_NEAR(geometric.modes[i].vacuum_wavenumber_squared/
            ref.vacuum_wavenumber_squared,0.25,1e-9);
        NOVA_REQUIRE_NEAR(magnetic.modes[i].vacuum_wavenumber_squared/
            ref.vacuum_wavenumber_squared,0.25,1e-9);
        NOVA_REQUIRE_NEAR(combined.modes[i].vacuum_wavenumber_squared/
            ref.vacuum_wavenumber_squared,0.0625,1e-9);
        NOVA_REQUIRE_NEAR(geometric.modes[i].frequency_hz/ref.frequency_hz,0.5,1e-9);
        NOVA_REQUIRE_NEAR(magnetic.modes[i].frequency_hz/ref.frequency_hz,0.5,1e-9);
        NOVA_REQUIRE_NEAR(combined.modes[i].frequency_hz/ref.frequency_hz,0.25,1e-9);
        NOVA_REQUIRE_NEAR(ref.frequency_hz/
            std::sqrt(ref.vacuum_wavenumber_squared),
            299792458.0/(2.0*std::numbers::pi),1e-7);
    }
}
