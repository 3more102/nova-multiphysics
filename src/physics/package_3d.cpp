#include "package_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace nova::physics {
namespace {
bool positive(double value) { return std::isfinite(value) && value>0; }

void validate(const PackageOptions3D& o) {
    if (!positive(o.minimum_dt) || !positive(o.initial_dt) || !positive(o.maximum_dt) ||
        o.minimum_dt>o.initial_dt || o.initial_dt>o.maximum_dt ||
        !positive(o.temperature_scale) || !positive(o.voltage_scale) ||
        !(o.error_kelvin>=0 && std::isfinite(o.error_kelvin)) ||
        !(o.relative_error>=0 && std::isfinite(o.relative_error)) ||
        !(o.relative_residual>=0 && o.relative_residual<1) ||
        !(o.relaxation>0 && o.relaxation<=1) || o.maximum_attempts==0 ||
        (o.strategy!=TransportStrategy3D::MonolithicNewton && o.strategy!=TransportStrategy3D::Segregated))
        throw std::invalid_argument("Invalid package time, scaling or nonlinear controls.");
    if (!std::is_sorted(o.events.begin(),o.events.end()) ||
        std::adjacent_find(o.events.begin(),o.events.end())!=o.events.end() ||
        !std::all_of(o.events.begin(),o.events.end(),[](double t){ return std::isfinite(t); }))
        throw std::invalid_argument("Package events must be finite and strictly increasing.");
}

void validate_state(const PackageState3D& s,std::size_t n) {
    if (!std::isfinite(s.time) || !std::isfinite(s.previous_time) ||
        !std::isfinite(s.next_dt) || s.next_dt<0 ||
        s.temperature.size()!=n || s.voltage.size()!=n ||
        !std::all_of(s.temperature.begin(),s.temperature.end(),positive) ||
        !std::all_of(s.voltage.begin(),s.voltage.end(),[](double v){ return std::isfinite(v); }) ||
        (!s.previous_temperature.empty() && (!(s.previous_time<s.time) ||
            s.previous_temperature.size()!=n ||
            !std::all_of(s.previous_temperature.begin(),s.previous_temperature.end(),positive))))
        throw std::invalid_argument("Invalid accepted package state or multistep history.");
}

std::vector<double> pack(const PackageState3D& s) {
    auto result=s.temperature;
    result.insert(result.end(),s.voltage.begin(),s.voltage.end());
    return result;
}

struct Trial {
    PackageState3D state;
    std::size_t iterations{};
    double residual{}, defect{}, reaction{};
    std::string failure;
};
}

PackageState3D PreparedPackage3D::initialize(std::vector<double> temperature,double t0) const {
    const auto n=model().mesh.nodes.size();
    if (temperature.size()!=n) throw std::invalid_argument("Initial package temperature size mismatch.");
    PackageState3D s;
    s.temperature=std::move(temperature); s.voltage.resize(n); s.time=t0; s.previous_time=t0;
    for (const auto& bc:model().temperature) s.temperature[bc.node]=bc.value(t0);
    for (const auto& bc:model().voltage) s.voltage[bc.node]=bc.value(t0);
    validate_state(s,n);
    impl_->electrical_predictor(s);
    return s;
}

void PreparedPackage3D::Impl::electrical_predictor(PackageState3D& s) const {
    const auto n=model.mesh.nodes.size();
    const auto sys=system(pack(s),s,s.time,{0,0,1},{},true);
    const auto diagonal=sys.jacobian.diagonal();
    std::vector<double> left(2*n,1), right(2*n,1), rhs(n);
    for (std::size_t i=0;i<n;++i) {
        left[n+i]=1/diagonal[n+i];
        rhs[i]=-sys.residual[n+i]*left[n+i];
    }
    linalg::GmresOptions linear;
    linear.solver={5000,1e-13,1e-15,linalg::PreconditionerType::IncompleteLU0};
    const auto solve=linalg::solve_gmres(sys.jacobian.scaled_block(n,n,left,right),rhs,linear);
    if (!solve.converged) throw std::runtime_error("Package initial electrical equilibrium: "+solve.message);
    for (std::size_t i=0;i<n;++i) s.voltage[i]+=solve.solution[i];
}

void PreparedPackage3D::advance(PackageState3D& state,double end_time,const PackageOptions3D& o,
    const std::function<void(const PackageState3D&,const PackageStep3D&)>& accepted) {
    const auto n=model().mesh.nodes.size();
    validate(o); validate_state(state,n);
    if (!std::isfinite(end_time) || end_time<state.time)
        throw std::invalid_argument("Package end time precedes the accepted state.");
    const auto step=[&](const PackageState3D& old,double t,time::ImplicitScheme scheme) {
        Trial trial;
        double previous_dt=old.previous_temperature.empty() ? 0 : old.time-old.previous_time;
        if (std::binary_search(o.events.begin(),o.events.end(),old.time)) previous_dt=0;
        const auto formula=time::implicit_formula(scheme,t-old.time,previous_dt);
        const auto past=impl_->old_operator(old,formula.theta);
        auto predicted=old;
        predicted.time=t;
        bool new_voltage=false;
        for (const auto& bc:model().voltage) {
            const double value=bc.value(t);
            new_voltage |= predicted.voltage[bc.node]!=value;
            predicted.voltage[bc.node]=value;
        }
        if (new_voltage) impl_->electrical_predictor(predicted);
        auto x=pack(predicted);
        for (const auto& bc:model().temperature) x[bc.node]=bc.value(t);
        const auto assemble=[&](const auto& v,bool tangent) {
            return impl_->system(v,old,t,formula,past,tangent);
        };
        auto initial=assemble(x,true);
        auto row=initial.jacobian.diagonal();
        std::vector<double> scale(2*n,o.voltage_scale);
        std::fill_n(scale.begin(),n,o.temperature_scale);
        double initial_norm=0;
        for (std::size_t i=0;i<2*n;++i) {
            const double magnitude=std::abs(row[i])*scale[i];
            if (!positive(magnitude)) { trial.failure="singular equation scale"; return trial; }
            row[i]=1/magnitude;
            initial_norm=std::hypot(initial_norm,row[i]*initial.residual[i]);
        }
        const double target=std::max(o.newton.residual_tolerance,o.relative_residual*initial_norm);
        bool first_tangent=true;
        // Solve increments about a physical anchor, so volts and kelvin have independent scales.
        const auto solve_block=[&](std::size_t first,std::size_t count,double relaxation) {
            const auto anchored=x;
            const auto physical=[&](std::span<const double> update) {
                auto v=anchored;
                for (std::size_t i=0;i<count;++i) v[first+i]+=scale[first+i]*update[i];
                return v;
            };
            const auto residual=[&](std::span<const double> u) {
                const auto v=physical(u);
                std::vector<double> r(count,std::numeric_limits<double>::infinity());
                if (!std::all_of(v.begin(),v.begin()+static_cast<std::ptrdiff_t>(n),positive)) return r;
                try {
                    const auto system=assemble(v,false);
                    for (std::size_t i=0;i<count;++i) r[i]=row[first+i]*system.residual[first+i];
                } catch (const std::domain_error&) { /* Armijo can shorten this material-domain excursion. */ }
                return r;
            };
            const auto jacobian=[&](std::span<const double> u) {
                if (first_tangent) {
                    first_tangent=false;
                    auto block=initial.jacobian.scaled_block(first,count,row,scale);
                    initial={};
                    return block;
                }
                return assemble(physical(u),true).jacobian.scaled_block(first,count,row,scale);
            };
            auto settings=o.newton;
            settings.relative_residual_tolerance=0;
            settings.residual_tolerance=count==n ? target*0.25 : target;
            const auto solve=nonlinear::solve_newton(residual,jacobian,std::vector<double>(count),settings);
            first_tangent=false;
            initial={};
            trial.iterations+=solve.iterations;
            if (!solve.converged) { trial.failure=solve.message; return; }
            for (std::size_t i=0;i<count;++i) x[first+i]+=relaxation*scale[first+i]*solve.solution[i];
        };
        auto check=[&]() {
            const auto sys=assemble(x,false);
            trial.residual=0;
            for (std::size_t i=0;i<2*n;++i) trial.residual=std::hypot(trial.residual,row[i]*sys.residual[i]);
            trial.defect=sys.energy_defect_watts; trial.reaction=sys.dirichlet_heat_watts;
            return trial.residual<=target;
        };
        if (o.strategy==TransportStrategy3D::MonolithicNewton) {
            solve_block(0,2*n,1);
            if (trial.failure.empty() && !check()) trial.failure="full residual exceeds tolerance";
        } else {
            for (std::size_t iter=0;iter<o.newton.max_iterations;++iter) {
                solve_block(n,n,1); // voltage at held temperature
                if (!trial.failure.empty()) break;
                solve_block(0,n,o.relaxation); // heat at held voltage
                if (!trial.failure.empty() || check()) break;
                if (iter+1==o.newton.max_iterations) trial.failure="segregated coupling iteration limit";
            }
        }
        if (trial.failure.empty()) {
            trial.state.time=t; trial.state.previous_time=old.time;
            trial.state.previous_temperature=old.temperature;
            trial.state.temperature.assign(x.begin(),x.begin()+static_cast<std::ptrdiff_t>(n));
            trial.state.voltage.assign(x.begin()+static_cast<std::ptrdiff_t>(n),x.end());
        }
        return trial;
    };

    const bool adaptive=o.error_kelvin>0 || o.relative_error>0;
    double dt=std::clamp(state.next_dt>0 ? state.next_dt : o.initial_dt,o.minimum_dt,o.maximum_dt);
    std::size_t attempts=0, rejected=0;
    while (state.time<end_time) {
        if (++attempts>o.maximum_attempts) throw std::runtime_error("Package attempt limit at last accepted state.");
        const auto event=std::upper_bound(o.events.begin(),o.events.end(),state.time);
        const double boundary=event==o.events.end() ? end_time : std::min(end_time,*event);
        if (o.scheme==time::ImplicitScheme::BDF2 && !state.previous_temperature.empty())
            dt=std::min(dt,2*(state.time-state.previous_time));
        // Snap only at roundoff relative to the time coordinates; t remains strictly increasing.
        const double ulp=8*std::numeric_limits<double>::epsilon()*std::max(std::abs(boundary),std::abs(state.time));
        const double t=(boundary-state.time<=dt+ulp) ? boundary : state.time+dt;
        if (!(t>state.time)) throw std::runtime_error("Package timestep cannot advance floating-point time.");
        const double h=t-state.time;
        auto high=step(state,t,o.scheme);
        double error=0;
        if (adaptive && high.failure.empty()) {
            Trial comparison;
            if (o.scheme==time::ImplicitScheme::BackwardEuler) {
                const double mid=state.time+h/2;
                if (!(mid>state.time && mid<t)) throw std::runtime_error("No representable half step.");
                auto half=step(state,mid,o.scheme);
                if (half.failure.empty()) {
                    comparison=step(half.state,t,o.scheme);
                    comparison.iterations+=half.iterations;
                    comparison.defect=(half.defect*(mid-state.time)+comparison.defect*(t-mid))/h;
                    comparison.reaction=(half.reaction*(mid-state.time)+comparison.reaction*(t-mid))/h;
                } else comparison=std::move(half);
                // Accept the refined BE trajectory and retain the outer accepted history.
                std::swap(high,comparison);
            } else comparison=step(state,t,time::ImplicitScheme::BackwardEuler);
            high.iterations+=comparison.iterations;
            if (!comparison.failure.empty()) high.failure=comparison.failure;
            if (high.failure.empty()) {
                high.state.previous_time=state.time;
                high.state.previous_temperature=state.temperature;
                for (std::size_t i=0;i<n;++i) {
                    const double T=high.state.temperature[i], other=comparison.state.temperature[i];
                    const double tol=o.error_kelvin+o.relative_error*std::max(std::abs(T),std::abs(other));
                    error=std::max(error,std::abs(T-other)/tol);
                }
            }
        }
        const double factor=adaptive && error>0 ? std::clamp(.9/std::sqrt(error),.2,2.) : 2.;
        PackageStep3D report;
        if (high.failure.empty() && error<=1) {
            report.mechanics=impl_->equilibrate(high.state.temperature,o.mechanics);
            if (!report.mechanics.linear_solver.converged) high.failure=report.mechanics.linear_solver.message;
        }
        if (!high.failure.empty() || error>1) {
            ++rejected;
            if (h<=o.minimum_dt) throw std::runtime_error("Package minimum step reached: "+
                (high.failure.empty() ? std::string("embedded temporal error") : high.failure));
            dt=std::max(o.minimum_dt,h*(high.failure.empty() ? std::min(.9,factor) : .5));
            continue;
        }
        report.transport=recover(high.state.temperature,high.state.voltage,t);
        report.residual=high.residual; report.normalized_error=error;
        report.energy_defect_watts=high.defect; report.dirichlet_heat_watts=high.reaction;
        report.rejected_before=rejected; report.iterations=high.iterations;
        dt=std::clamp(h*factor,o.minimum_dt,adaptive ? o.maximum_dt : o.initial_dt);
        high.state.next_dt=dt;
        state=std::move(high.state);
        if (accepted) accepted(state,report);
        rejected=0;
    }
}
}
