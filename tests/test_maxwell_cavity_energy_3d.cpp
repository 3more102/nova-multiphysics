#include "test_framework.hpp"
#include "nova/physics/maxwell_cavity_energy_3d.hpp"
#include <algorithm>
#include <stdexcept>

using namespace nova::physics;

NOVA_TEST(cavity_energy_single_tetra_exact_affine_integrals) {
    const nova::fem::TetrahedralMesh3D mesh{
        {{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}}};
    const std::vector<MaxwellCellMaterial3D> materials(1);
    MaxwellCavityMatrices3D matrices;
    const auto tet=nova::fem::make_nedelec_tetrahedron_3d(mesh,0);
    matrices.edges.assign(tet.edges.begin(),tet.edges.end());
    matrices.free_edges={0,1,2,3,4,5};
    MaxwellCavityMode3D mode;
    mode.vacuum_wavenumber_squared=16.0;
    mode.edge_field={1,0,0,0,0,0};
    const auto r=evaluate_maxwell_cavity_energy_3d(mesh,materials,matrices,mode);
    NOVA_REQUIRE_NEAR(r.electric_integral,1.0/12.0,1e-14);
    NOVA_REQUIRE_NEAR(r.magnetic_integral,1.0/12.0,1e-14);
    NOVA_REQUIRE_NEAR(r.peak_electric_weight,3.0,1e-14);
    NOVA_REQUIRE_NEAR(r.electric_mode_volume_m3,1.0/36.0,1e-14);
    NOVA_REQUIRE_NEAR(r.cells[0].total_participation,1.0,1e-14);
}

NOVA_TEST(cavity_energy_physical_mode_equipartition_and_amplitude_invariance) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
    for (std::size_t i=0;i<materials.size();++i)
        materials[i].relative_permittivity=i%2==0?2.0:3.0;
    const auto matrices=assemble_maxwell_cavity_3d(mesh,materials);
    const auto solved=solve_maxwell_cavity_3d(matrices,2);
    NOVA_REQUIRE(solved.converged);
    NOVA_REQUIRE(solved.modes.size()==2);
    for (auto mode:solved.modes) {
        const auto r=evaluate_maxwell_cavity_energy_3d(mesh,materials,matrices,mode);
        NOVA_REQUIRE_NEAR(r.electric_integral,1.0,1e-10);
        NOVA_REQUIRE(r.relative_energy_imbalance<1e-9);
        NOVA_REQUIRE(r.electric_mode_volume_m3>0.0 && r.electric_mode_volume_m3<=1.0);
        double electric=0,total=0;
        for (const auto& c:r.cells) {
            NOVA_REQUIRE(c.electric_participation>=0.0 && c.total_participation>=0.0);
            electric+=c.electric_participation;
            total+=c.total_participation;
        }
        NOVA_REQUIRE_NEAR(electric,1.0,1e-12);
        NOVA_REQUIRE_NEAR(total,1.0,1e-12);
        for (auto& e:mode.edge_field) e*=-3;
        const auto scaled=evaluate_maxwell_cavity_energy_3d(mesh,materials,matrices,mode);
        NOVA_REQUIRE_NEAR(scaled.electric_integral,9*r.electric_integral,1e-10);
        NOVA_REQUIRE_NEAR(scaled.electric_mode_volume_m3,r.electric_mode_volume_m3,1e-12);
        for (std::size_t i=0;i<r.cells.size();++i)
            NOVA_REQUIRE_NEAR(scaled.cells[i].total_participation,r.cells[i].total_participation,1e-12);
    }
}

NOVA_TEST(cavity_energy_rejects_invalid_mode_and_material) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
    const auto matrices=assemble_maxwell_cavity_3d(mesh,materials);
    const auto solved=solve_maxwell_cavity_3d(matrices,1);
    NOVA_REQUIRE(solved.converged && solved.modes.size()==1);
    auto mode=solved.modes[0];
    const auto rejects=[&](){
        try {(void)evaluate_maxwell_cavity_energy_3d(mesh,materials,matrices,mode);}
        catch (const std::invalid_argument&) {return true;}
        return false;
    };
    materials[0].relative_permittivity={1.0,0.1};
    NOVA_REQUIRE(rejects());
    materials[0].relative_permittivity=1.0;
    mode.vacuum_wavenumber_squared=0;
    NOVA_REQUIRE(rejects());
    mode=solved.modes[0];
    std::fill(mode.edge_field.begin(),mode.edge_field.end(),0.0);
    NOVA_REQUIRE(rejects());
}
