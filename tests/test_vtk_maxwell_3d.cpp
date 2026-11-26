#include "test_framework.hpp"
#include "nova/io/vtk_maxwell_3d.hpp"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace nova::fem;
using namespace nova::physics;
using namespace nova::io;

NOVA_TEST(maxwell_cell_output_affine_field_and_exact_loss) {
    using C=std::complex<double>;
    TetrahedralMesh3D mesh{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{2,0,3,1}}};
    const std::array<MaxwellCellMaterial3D,1> materials{{{{2,0.4},1}}};
    const auto sys=assemble_maxwell_3d(mesh,materials,2,false);
    std::vector<C> edge;
    // E=(-y,x,i), exact circulation on each globally oriented edge.
    for (const auto& e : sys.edges) {
        const auto a=mesh.nodes[e[0]], b=mesh.nodes[e[1]];
        edge.push_back(-(a.y+b.y)/2*(b.x-a.x)+(a.x+b.x)/2*(b.y-a.y)+C(0,b.z-a.z));
    }
    const auto fields=reconstruct_maxwell_cells_3d(mesh,materials,2,sys,edge);
    NOVA_REQUIRE(fields.size()==1);
    const auto& f=fields.front();
    NOVA_REQUIRE_NEAR(std::abs(f.centroid.electric_field[0]+0.25),0,1e-14);
    NOVA_REQUIRE_NEAR(std::abs(f.centroid.electric_field[1]-0.25),0,1e-14);
    NOVA_REQUIRE_NEAR(std::abs(f.centroid.electric_field[2]-C(0,1)),0,1e-14);
    constexpr double factor=0.5*2*299792458.0*8.8541878128e-12*0.4;
    NOVA_REQUIRE_NEAR(f.mean_absorbed_power_density,factor*1.2,1e-14);
    NOVA_REQUIRE_NEAR(f.absorbed_power,factor*1.2/6,1e-14);
    NOVA_REQUIRE(f.mean_absorbed_power_density>f.centroid.absorbed_power_density);
    std::ostringstream vtk;
    vtk<<std::fixed<<std::setprecision(3);
    const auto flags=vtk.flags();
    write_legacy_vtk_maxwell_3d(vtk,mesh,fields);
    NOVA_REQUIRE(vtk.flags()==flags && vtk.precision()==3);
    const auto content=vtk.str();
    NOVA_REQUIRE(content.find("CELL_DATA 1\n")!=std::string::npos);
    NOVA_REQUIRE(content.find("POINT_DATA")==std::string::npos);
    const std::string key="VECTORS electric_real double\n";
    const auto pos=content.find(key);
    NOVA_REQUIRE(pos!=std::string::npos);
    std::istringstream data(content.substr(pos+key.size()));
    double x{},y{},z{};
    NOVA_REQUIRE(static_cast<bool>(data>>x>>y>>z));
    NOVA_REQUIRE_NEAR(x,-0.25,1e-14);
    NOVA_REQUIRE_NEAR(y,0.25,1e-14);
    NOVA_REQUIRE_NEAR(z,0,1e-14);
    NOVA_REQUIRE(content.find("VECTORS electric_imag double")!=std::string::npos);
    NOVA_REQUIRE(content.find("SCALARS mean_absorbed_power_density double 1")!=std::string::npos);
}

NOVA_TEST(maxwell_cell_output_multicell_and_invalid_inputs) {
    const auto mesh=make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size(),{{2,0},1});
    auto sys=assemble_maxwell_3d(mesh,materials,2,false);
    std::vector<std::complex<double>> edge;
    for (const auto& e:sys.edges)
        edge.emplace_back(mesh.nodes[e[1]].x-mesh.nodes[e[0]].x,0);
    auto fields=reconstruct_maxwell_cells_3d(mesh,materials,2,sys,edge);
    NOVA_REQUIRE(fields.size()==6);
    for (const auto& f:fields) {
        NOVA_REQUIRE_NEAR(f.centroid.electric_field[0].real(),1,1e-14);
        NOVA_REQUIRE_NEAR(f.absorbed_power,0,0);
    }
    std::swap(sys.edges[0],sys.edges[1]);
    bool rejected=false;
    try { (void)reconstruct_maxwell_cells_3d(mesh,materials,2,sys,edge); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
    std::ostringstream vtk;
    fields.back().mean_absorbed_power_density=std::numeric_limits<double>::infinity();
    rejected=false;
    try { write_legacy_vtk_maxwell_3d(vtk,mesh,fields); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected && vtk.str().empty());
    fields.back().mean_absorbed_power_density=0;
    vtk.setstate(std::ios::badbit);
    rejected=false;
    try { write_legacy_vtk_maxwell_3d(vtk,mesh,fields); }
    catch (const std::runtime_error&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}
