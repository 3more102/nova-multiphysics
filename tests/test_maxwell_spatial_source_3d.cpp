#include "test_framework.hpp"
#include "nova/physics/maxwell_spatial_source_3d.hpp"

#include <limits>

using namespace nova::physics;
using namespace nova::fem;

namespace {
TetrahedralMesh3D source_unit_cell() {
    return {{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}}};
}
}

NOVA_TEST(maxwell_spatial_source_matches_heterogeneous_constant_load) {
    const auto mesh=make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size(),{{2.0,0.4},1.0});
    const auto sys=assemble_maxwell_3d(mesh,materials,2.0,true);
    std::vector<ComplexVector3D> sources(mesh.elements.size());
    for (std::size_t c=0;c<sources.size();++c)
        sources[c]={std::complex<double>(0.1*c,0.3),{2.0,-0.5},{-0.2,0.0}};
    const auto exact=assemble_maxwell_cell_load_3d(mesh,sources,sys);
    const auto actual=assemble_maxwell_spatial_load_3d(mesh,
        [&](std::size_t c,const auto&){return sources[c];},sys);
    NOVA_REQUIRE(actual.size()==exact.size());
    for (std::size_t e=0;e<exact.size();++e)
        NOVA_REQUIRE_NEAR(std::abs(actual[e]-exact[e]),0.0,1e-12);
}

NOVA_TEST(maxwell_spatial_source_quartic_exact_and_orientation_invariant) {
    auto mesh=source_unit_cell();
    const std::array<MaxwellCellMaterial3D,1> materials{{{{2.0,0.4},1.0}}};
    const auto sys=assemble_maxwell_3d(mesh,materials,2.0,false);
    const auto source=[](std::size_t,const std::array<double,3>& p) {
        const double x2=p[0]*p[0];
        return ComplexVector3D{std::complex<double>(1.0,0.5)*x2*x2,0.0,0.0};
    };
    const auto a=assemble_maxwell_spatial_load_3d(mesh,source,sys);
    // Edge 01: basis x-component 1-y-z. Edges 02,03: y,z.
    // Edges 12,13: -y,-z. Edge 23: zero.
    const std::array<double,6> exact{1.0/280,1.0/1680,1.0/1680,
        -1.0/1680,-1.0/1680,0.0};
    NOVA_REQUIRE(a.size()==6);
    for (std::size_t e=0;e<6;++e)
        NOVA_REQUIRE_NEAR(std::abs(a[e]-std::complex<double>(1.0,0.5)*exact[e]),0.0,1e-14);
    mesh.elements[0]={2,0,3,1};
    const auto b=assemble_maxwell_spatial_load_3d(mesh,source,sys);
    for (std::size_t e=0;e<6;++e)
        NOVA_REQUIRE_NEAR(std::abs(a[e]-b[e]),0.0,1e-14);
}

NOVA_TEST(maxwell_spatial_source_invalid_inputs) {
    const auto mesh=source_unit_cell();
    const std::array<MaxwellCellMaterial3D,1> materials{{{{2.0,0.4},1.0}}};
    auto sys=assemble_maxwell_3d(mesh,materials,2.0,false);
    const auto good=[](std::size_t,const auto&){return ComplexVector3D{};};
    const auto bad=[](std::size_t,const auto&){
        return ComplexVector3D{std::numeric_limits<double>::infinity(),0.0,0.0};
    };
    bool rejected=false;
    try {(void)assemble_maxwell_spatial_load_3d(mesh,{},sys);}
    catch(const std::invalid_argument&) {rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    try {(void)assemble_maxwell_spatial_load_3d(mesh,bad,sys);}
    catch(const std::invalid_argument&) {rejected=true;}
    NOVA_REQUIRE(rejected);
    std::swap(sys.edges[0],sys.edges[1]);
    rejected=false;
    try {(void)assemble_maxwell_spatial_load_3d(mesh,good,sys);}
    catch(const std::invalid_argument&) {rejected=true;}
    NOVA_REQUIRE(rejected);
}
