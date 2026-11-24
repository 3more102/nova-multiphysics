#pragma once
#include "nova/physics/package_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

namespace nova::physics {
struct PreparedPackage3D::Impl {
    struct Cell {
        std::array<ElasticVector3D,4> gradient;
        double volume{};
    };
    PackageModel3D model;
    std::vector<Cell> cells;
    linalg::CsrAssemblyPattern pattern;
    std::vector<double> areas;
    std::vector<unsigned char> active, fixed;
    std::unique_ptr<PreparedElasticity3D> elasticity;
    bool variable_modulus{};
    std::vector<ThermalExpansion3D> expansion;

    explicit Impl(PackageModel3D input);
    struct Evaluation {
        std::vector<double> residual;
        std::vector<double> entries;
    };
    void integrate(std::span<const double> temperature, std::span<const double> voltage, double t,
        const PackageState3D* old, time::ImplicitFormula formula, bool tangent,
        bool electric_rows, Evaluation& out, TransportFields3D* output = nullptr) const;
    TransportAssembly3D system(std::span<const double> x, const PackageState3D& old,
        double t, time::ImplicitFormula formula, std::span<const double> old_spatial, bool tangent) const;
    std::vector<double> old_operator(const PackageState3D& old, double theta) const;
    void electrical_predictor(PackageState3D& state) const;
    LinearElasticity3DResult equilibrate(std::span<const double> temperature,
        const linalg::SolverOptions& options);
};
}
