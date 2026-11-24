#include "test_framework.hpp"
#include "nova/physics/maxwell_3d.hpp"

#include <algorithm>
#include <complex>
#include <limits>
#include <stdexcept>

using namespace nova::fem;
using namespace nova::physics;

namespace {
TetrahedralMesh3D unit_cell() {
    return {{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}}};
}
}

NOVA_TEST(nedelec_exact_unit_integrals_and_gradient_kernel) {
    const auto cell=make_nedelec_tetrahedron_3d(unit_cell(),0);
    const auto m=nedelec_element_matrices_3d(cell);
    NOVA_REQUIRE_NEAR(cell.volume,1.0/6.0,1e-15);
    NOVA_REQUIRE_NEAR(m.mass[0][0],1.0/12.0,1e-15);
    NOVA_REQUIRE_NEAR(m.curl_curl[0][0],4.0/3.0,1e-15);
    const std::array<double,4> potential{0.2,0.7,-0.8,1.1};
    std::array<double,6> grad{};
    for (std::size_t e=0;e<6;++e)
        grad[e]=potential[cell.edges[e][1]]-potential[cell.edges[e][0]];
    for (std::size_t e=0;e<6;++e) {
        double v=0;
        for (std::size_t f=0;f<6;++f) v+=m.curl_curl[e][f]*grad[f];
        NOVA_REQUIRE_NEAR(v,0.0,2e-15);
    }
}

NOVA_TEST(nedelec_oriented_edge_circulation) {
    auto mesh=unit_cell();
    mesh.elements[0]={2,0,3,1};
    const auto cell=make_nedelec_tetrahedron_3d(mesh,0);
    for (std::size_t e=0;e<6;++e) {
        const auto [i,j]=tetrahedron_local_edges[e];
        std::array<double,4> l{};
        l[i]=l[j]=0.5;
        const auto b=evaluate_nedelec_basis_3d(cell,l[1],l[2],l[3]);
        const auto a=mesh.nodes[cell.edges[e][0]], z=mesh.nodes[cell.edges[e][1]];
        for (std::size_t f=0;f<6;++f) {
            const auto& w=b[f];
            NOVA_REQUIRE_NEAR(w[0]*(z.x-a.x)+w[1]*(z.y-a.y)+w[2]*(z.z-a.z),
                              e==f ? 1.0 : 0.0,1e-14);
        }
    }
}

NOVA_TEST(maxwell_permutation_invariant_complex_assembly) {
    auto mesh=unit_cell();
    const std::array<MaxwellCellMaterial3D,1> mat{{{{2.0,0.3},1.2}}};
    const auto a=assemble_maxwell_3d(mesh,mat,1.7,false);
    mesh.elements[0]={2,0,3,1};
    const auto b=assemble_maxwell_3d(mesh,mat,1.7,false);
    NOVA_REQUIRE(a.edges==b.edges);
    for (std::size_t i=0;i<a.real_block.rows();++i)
        for (std::size_t j=0;j<a.real_block.columns();++j)
            NOVA_REQUIRE_NEAR(a.real_block.at(i,j),b.real_block.at(i,j),1e-14);
}

NOVA_TEST(maxwell_lossy_constant_field_manufactured_solution) {
    using C=std::complex<double>;
    const ComplexVector3D exact{C(1.0,0.2),C(-0.3,0.4),C(0.5,-0.1)};
    for (std::size_t n : {1U,2U}) {
        const auto mesh=make_box_mesh_3d(0,1,0,1,0,1,n,n,n).mesh;
        const C eps(2.0,0.6);
        const double k=2.3;
        const std::vector<MaxwellCellMaterial3D> mat(mesh.elements.size(),{eps,1.0});
        ComplexVector3D f{};
        for (std::size_t d=0;d<3;++d) f[d]=-k*k*eps*exact[d];
        const std::vector<ComplexVector3D> sources(mesh.elements.size(),f);
        const auto sys=assemble_maxwell_3d(mesh,mat,k,false);
        const auto load=assemble_maxwell_cell_load_3d(mesh,sources,sys);
        nova::linalg::GmresOptions opts;
        opts.restart=100;
        opts.solver.max_iterations=1000;
        opts.solver.relative_tolerance=1e-12;
        const auto result=solve_maxwell_3d(sys,load,opts);
        NOVA_REQUIRE(result.diagnostics.converged);
        double max_error=0.0;
        for (std::size_t e=0;e<sys.edges.size();++e) {
            const auto a=mesh.nodes[sys.edges[e][0]], b=mesh.nodes[sys.edges[e][1]];
            const C expected=exact[0]*(b.x-a.x)+exact[1]*(b.y-a.y)+exact[2]*(b.z-a.z);
            max_error=std::max(max_error,std::abs(result.edge_field[e]-expected));
        }
        std::cout<<"Maxwell constant field: cells="<<mesh.elements.size()
                 <<" max circulation error="<<max_error
                 <<" residual="<<result.diagnostics.residual_norm<<'\n';
        NOVA_REQUIRE_NEAR(max_error,0.0,1e-9);
    }
}

NOVA_TEST(maxwell_pec_boundary_elimination_and_zero_solution) {
    const auto mesh=make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    const std::vector<MaxwellCellMaterial3D> mat(mesh.elements.size(),{{2.0,0.5},1.0});
    const auto sys=assemble_maxwell_3d(mesh,mat,2.0,true);
    NOVA_REQUIRE(sys.free_edges.size()==1);
    NOVA_REQUIRE(sys.real_block.rows()==2);
    const auto result=solve_maxwell_3d(sys,std::vector<std::complex<double>>(sys.edges.size()));
    NOVA_REQUIRE(result.diagnostics.converged);
    for (auto z : result.edge_field) NOVA_REQUIRE_NEAR(std::abs(z),0.0,1e-15);
}

NOVA_TEST(maxwell_rejects_active_material_and_invalid_frequency) {
    const auto mesh=unit_cell();
    std::array<MaxwellCellMaterial3D,1> mat{{{{2.0,-0.1},1.0}}};
    bool rejected=false;
    try { (void)assemble_maxwell_3d(mesh,mat,2.0,false); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
    mat[0].relative_permittivity={2.0,0.1};
    rejected=false;
    try { (void)assemble_maxwell_3d(mesh,mat,std::numeric_limits<double>::infinity(),false); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}

NOVA_TEST(maxwell_affine_field_curl_poynting_and_absorption) {
    using C=std::complex<double>;
    const auto mesh=unit_cell();
    const auto cell=make_nedelec_tetrahedron_3d(mesh,0);
    const MaxwellCellMaterial3D material{{2.0,0.4},1.0};
    const double k=2.0;
    constexpr double c=299792458.0, eps0=8.8541878128e-12;
    const double omega=k*c;
    const double mu0=1.0/(eps0*c*c);
    // E=(-y,x,i), curl(E)=(0,0,2); integral along each edge is
    // exactly the midpoint value dotted with the oriented edge vector.
    std::array<C,6> coefficients{};
    for (std::size_t e=0;e<6;++e) {
        const auto a=mesh.nodes[cell.edges[e][0]], b=mesh.nodes[cell.edges[e][1]];
        coefficients[e]=-(a.y+b.y)/2.0*(b.x-a.x)
                       +(a.x+b.x)/2.0*(b.y-a.y)+C(0,b.z-a.z);
    }
    const auto s=evaluate_maxwell_field_3d(cell,coefficients,material,k,0.2,0.3,0.1);
    NOVA_REQUIRE_NEAR(std::abs(s.electric_field[0]+0.3),0.0,1e-14);
    NOVA_REQUIRE_NEAR(std::abs(s.electric_field[1]-0.2),0.0,1e-14);
    NOVA_REQUIRE_NEAR(std::abs(s.electric_field[2]-C(0,1)),0.0,1e-14);
    NOVA_REQUIRE_NEAR(std::abs(s.magnetic_field[2]-C(0,-2/(omega*mu0))),0.0,1e-14);
    NOVA_REQUIRE_NEAR(s.time_averaged_poynting[0],0.0,1e-14);
    NOVA_REQUIRE_NEAR(s.time_averaged_poynting[1],0.0,1e-14);
    NOVA_REQUIRE_NEAR(s.time_averaged_poynting[2],0.0,1e-14);
    NOVA_REQUIRE_NEAR(s.absorbed_power_density,0.5*omega*eps0*0.4*1.13,1e-14);
    const std::array<MaxwellCellMaterial3D,1> materials{material};
    const auto sys=assemble_maxwell_3d(mesh,materials,k,false);
    NOVA_REQUIRE(std::equal(sys.edges.begin(),sys.edges.end(),cell.edges.begin(),cell.edges.end()));
    const auto q=integrate_maxwell_absorption_3d(mesh,materials,k,sys,coefficients);
    // Integral of x^2+y^2+1 over unit tetrahedron is V*(1+2/10).
    const double exact=0.5*omega*eps0*0.4*1.2/6.0;
    NOVA_REQUIRE_NEAR(q[0],exact,1e-14);
    std::cout<<"Maxwell dielectric absorption: numerical="<<q[0]
             <<" exact="<<exact<<" absolute error="<<std::abs(q[0]-exact)<<'\n';
}

NOVA_TEST(maxwell_poynting_nonzero_and_lossless_material) {
    using C=std::complex<double>;
    const auto mesh=unit_cell();
    const auto cell=make_nedelec_tetrahedron_3d(mesh,0);
    // E=(1-i*y,i*x,0), curl(E)=(0,0,2i), so H_z is positive real.
    std::array<C,6> coefficients{};
    for (std::size_t e=0;e<6;++e) {
        const auto a=mesh.nodes[cell.edges[e][0]], b=mesh.nodes[cell.edges[e][1]];
        coefficients[e]=C(1,-(a.y+b.y)/2.0)*(b.x-a.x)
                       +C(0,(a.x+b.x)/2.0)*(b.y-a.y);
    }
    const auto s=evaluate_maxwell_field_3d(cell,coefficients,{{2,0},1},2,0.2,0.3,0.1);
    NOVA_REQUIRE(s.magnetic_field[2].real()>0.0);
    NOVA_REQUIRE_NEAR(s.time_averaged_poynting[1],-0.5*s.magnetic_field[2].real(),1e-14);
    NOVA_REQUIRE_NEAR(s.absorbed_power_density,0.0,0.0);
}
