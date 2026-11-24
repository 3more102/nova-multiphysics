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

NOVA_TEST(maxwell_uniaxial_cavity_exact_one_dof_mass_and_frequency) {
    using namespace nova::fem;
    const auto mesh=make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    const double invsqrt2=1.0/std::sqrt(2.0);
    const MaxwellUniaxialMaterial3D material{{2,0},{5,0},{1,1,0},1.4};
    const std::vector<MaxwellUniaxialMaterial3D> materials(mesh.elements.size(),material);
    const auto a=assemble_maxwell_uniaxial_cavity_3d(mesh,materials);
    NOVA_REQUIRE(a.free_edges.size()==1);
    const auto edge=a.edges[a.free_edges[0]];
    constexpr std::array<std::array<double,3>,4> vertices{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,1}
    }};
    double expected_mass=0.0,expected_stiffness=0.0;
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=make_nedelec_tetrahedron_3d(mesh,c);
        const auto it=std::find(cell.edges.begin(),cell.edges.end(),edge);
        if (it==cell.edges.end()) continue;
        const auto e=static_cast<std::size_t>(it-cell.edges.begin());
        std::array<EdgeVector3D,4> v{};
        for (std::size_t i=0;i<4;++i)
            v[i]=evaluate_nedelec_basis_3d(cell,vertices[i][0],vertices[i][1],vertices[i][2])[e];
        // Integrate affine products analytically with barycentric moments,
        // not the four-point rule used by the assembler.
        for (std::size_t i=0;i<4;++i)
            for (std::size_t j=0;j<4;++j) {
                double dot=0.0;
                for (std::size_t d=0;d<3;++d) dot+=v[i][d]*v[j][d];
                const double pi=(v[i][0]+v[i][1])*invsqrt2;
                const double pj=(v[j][0]+v[j][1])*invsqrt2;
                expected_mass+=cell.volume*(i==j?2.0:1.0)/20.0*(2.0*dot+3.0*pi*pj);
            }
        for (double curl : cell.curls[e])
            expected_stiffness+=cell.volume*curl*curl/1.4;
    }
    NOVA_REQUIRE_NEAR(a.mass.at(0,0),expected_mass,1e-13);
    NOVA_REQUIRE_NEAR(a.stiffness.at(0,0),expected_stiffness,1e-13);
    const auto result=solve_maxwell_cavity_3d(a,1);
    NOVA_REQUIRE(result.converged && result.modes.size()==1);
    NOVA_REQUIRE(result.nullspace_dimension==0);
    NOVA_REQUIRE_NEAR(result.modes[0].vacuum_wavenumber_squared,
        expected_stiffness/expected_mass,1e-12);
}

NOVA_TEST(maxwell_uniaxial_cavity_isotropic_limit_and_spectral_bounds) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    std::vector<MaxwellUniaxialMaterial3D> mat(mesh.elements.size(),{{2,0},{2,0},{1,2,3},1});
    const std::vector<MaxwellCellMaterial3D> iso(mesh.elements.size(),{{2,0},1});
    const auto a=assemble_maxwell_cavity_3d(mesh,iso);
    const auto b=assemble_maxwell_uniaxial_cavity_3d(mesh,mat);
    NOVA_REQUIRE(a.edges==b.edges && a.free_edges==b.free_edges);
    for (std::size_t i=0;i<a.mass.rows();++i)
        for (std::size_t j=0;j<a.mass.columns();++j) {
            NOVA_REQUIRE_NEAR(a.mass.at(i,j),b.mass.at(i,j),1e-13);
            NOVA_REQUIRE_NEAR(a.stiffness.at(i,j),b.stiffness.at(i,j),1e-13);
        }
    const auto baseline=solve_maxwell_cavity_3d(a,3);
    for (auto& m : mat) m.extraordinary_permittivity={6,0};
    const auto tensor=assemble_maxwell_uniaxial_cavity_3d(mesh,mat);
    const auto modes=solve_maxwell_cavity_3d(tensor,3);
    NOVA_REQUIRE(baseline.converged && modes.converged);
    NOVA_REQUIRE(baseline.modes.size()==3 && modes.modes.size()==3);
    NOVA_REQUIRE(modes.nullspace_dimension==baseline.nullspace_dimension);
    // 2 I <= epsilon <= 6 I gives baseline/3 <= lambda <= baseline.
    for (std::size_t i=0;i<3;++i) {
        const double ref=baseline.modes[i].vacuum_wavenumber_squared;
        const double lambda=modes.modes[i].vacuum_wavenumber_squared;
        NOVA_REQUIRE(lambda>=ref/3.0-1e-8 && lambda<=ref+1e-8);
        NOVA_REQUIRE(modes.modes[i].relative_residual<1e-10);
    }
    const double k=2.3;
    const auto driven=assemble_maxwell_uniaxial_3d(mesh,mat,k,true);
    for (std::size_t i=0;i<tensor.mass.rows();++i)
        for (std::size_t j=0;j<tensor.mass.columns();++j)
            NOVA_REQUIRE_NEAR(driven.real_block.at(i,j),
                tensor.stiffness.at(i,j)-k*k*tensor.mass.at(i,j),1e-13);
}

NOVA_TEST(maxwell_uniaxial_cavity_rejects_loss_and_indefinite_mass) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    std::vector<MaxwellUniaxialMaterial3D> mat(mesh.elements.size());
    for (const auto eps : {std::complex<double>{2,0.1},std::complex<double>{-1,0},
                           std::complex<double>{0,0}}) {
        mat[0].extraordinary_permittivity=eps;
        bool rejected=false;
        try {(void)assemble_maxwell_uniaxial_cavity_3d(mesh,mat);}
        catch(const std::invalid_argument&){rejected=true;}
        NOVA_REQUIRE(rejected);
    }
}
