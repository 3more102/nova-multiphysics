#include "nova/physics/optical_helmholtz_power_audit.hpp"
#include <iostream>
#include <stdexcept>

using namespace nova::optics;
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}
int main() {
    try {
        using C = std::complex<double>;
        constexpr double wavelength = 1.55e-6;
        std::vector<double> nodes{0.0};
        std::vector<C> eps;
        for (int j=0; j<180; ++j) {
            nodes.push_back(nodes.back()+wavelength*(j%2==0 ? 0.003 : 0.005));
            eps.push_back(j<60 ? C{2.25,0} : (j<120 ? C{3.2,0.8} : C{-4,1.2}));
        }
        const auto solution=solve_helmholtz_1d(wavelength,1.2,1.5,nodes,eps);
        const auto audit=audit_helmholtz_cell_power(wavelength,1.2,nodes,eps,solution);
        require(audit.cells.size()==eps.size(),"audit cell count");
        require(audit.max_interface_jump<1e-9,"interface continuity");
        require(std::abs(audit.entrance_residual)<1e-9,"incident minus reflected power");
        require(std::abs(audit.exit_residual)<1e-9,"transmitted power");
        for (std::size_t j=0; j<eps.size(); ++j) {
            const auto& c=audit.cells[j];
            require(std::abs(c.absorption_residual)<1e-9,"local absorption identity");
            require(c.left_flux>=c.right_flux-1e-9,"passive local power");
        }
        require(audit.cells[60].left_flux>audit.cells[119].right_flux,"lossy region dissipates");
        auto corrupt=solution;
        corrupt.electric_field[100]*=1.1;
        require(audit_helmholtz_cell_power(wavelength,1.2,nodes,eps,corrupt).max_interface_jump>1e-6,
                "detect corrupted field");
        corrupt=solution;
        corrupt.cell_absorptance[100]+=0.01;
        require(std::abs(audit_helmholtz_cell_power(wavelength,1.2,nodes,eps,corrupt).cells[100].absorption_residual)>0.009,
                "detect corrupted loss");
        bool rejected=false;
        try { (void)audit_helmholtz_cell_power(wavelength,0.0,nodes,eps,solution); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject invalid incident index");
        std::cout<<"NOVA local Helmholtz power audit tests passed\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
