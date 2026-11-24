#include "nova/io/package_study_3d.hpp"
#include "nova/io/thermal_checkpoint_3d.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace nova::io {
namespace {
void validate(const physics::PackageState3D& s,std::size_t n,const std::string& id) {
    const auto kelvin=[](double t){return std::isfinite(t) && t>0;};
    if (id.empty() || id.size()>128 || id.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-")!=std::string::npos ||
        !std::isfinite(s.time) || !std::isfinite(s.previous_time) || !std::isfinite(s.next_dt) || s.next_dt<=0 ||
        s.temperature.size()!=n || s.voltage.size()!=n || !std::all_of(s.temperature.begin(),s.temperature.end(),kelvin) ||
        !std::all_of(s.voltage.begin(),s.voltage.end(),[](double v){return std::isfinite(v);}) ||
        s.previous_temperature.size()!=n || !(s.previous_time<s.time) ||
        !std::all_of(s.previous_temperature.begin(),s.previous_temperature.end(),kelvin))
        throw std::invalid_argument("Invalid accepted package checkpoint state or model identifier.");
}
}
void write_package_checkpoint_3d(std::ostream& out,const fem::TetrahedralMesh3D& mesh,
    const std::string& id,const physics::PackageState3D& s) {
    validate(s,mesh.nodes.size(),id);
    std::ostringstream record;
    record.imbue(std::locale::classic());
    record<<std::setprecision(std::numeric_limits<double>::max_digits10)
          <<"NOVA_PACKAGE_CHECKPOINT 1 "<<id<<'\n'<<s.previous_time<<' '<<s.next_dt<<' '<<s.voltage.size()<<'\n';
    for (std::size_t i=0;i<s.voltage.size();++i) record<<s.voltage[i]<<' '<<s.previous_temperature[i]<<'\n';
    write_transient_heat_checkpoint_3d(record,mesh,s.time,s.temperature);
    out<<record.str();
    if (!out) throw std::runtime_error("Failed to write package checkpoint.");
}
physics::PackageState3D read_package_checkpoint_3d(std::istream& in,const fem::TetrahedralMesh3D& mesh,const std::string& id) {
    std::string magic,version,stored;
    physics::PackageState3D s;
    std::size_t n{};
    if (!(in>>magic>>version>>stored) || magic!="NOVA_PACKAGE_CHECKPOINT" || version!="1" || stored!=id ||
        !(in>>s.previous_time>>s.next_dt>>n) || n!=mesh.nodes.size())
        throw std::invalid_argument("Package checkpoint header, mesh size or model identity mismatch.");
    s.voltage.resize(n); s.previous_temperature.resize(n);
    for (std::size_t i=0;i<n;++i) if (!(in>>s.voltage[i]>>s.previous_temperature[i]))
        throw std::invalid_argument("Incomplete package checkpoint history.");
    auto thermal=read_transient_heat_checkpoint_3d(in,mesh);
    s.time=thermal.time; s.temperature=std::move(thermal.temperature);
    validate(s,n,id);
    return s;
}
}
