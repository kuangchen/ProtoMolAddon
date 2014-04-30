#include <fstream>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/addon/ion_trap/LQT.h>
#include <protomol/addon/util/SIAtomProxy.h>

namespace pt = boost::property_tree;
using namespace ProtoMolAddon::IonTrap;

LQT::LQTSpec::LQTSpec(const string &fname) {
  pt::ptree tree;
  pt::read_xml(fname, tree);

  double r0 = tree.get<double>("LQTSpec.r0");
  double z0 = tree.get<double>("LQTSpec.z0");
  double v_rf = tree.get<double>("LQTSpec.v_rf");
  double v_ec = tree.get<double>("LQTSpec.v_ec");
  double kappa = tree.get<double>("LQTSpec.kappa");
  double omega = tree.get<double>("LQTSpec.omega");
  (*this) = LQTSpec(r0, z0, omega, v_rf, v_ec, kappa);
}

LQT::LQT(const LQT::LQTSpec &spec) :
  spec(spec), 
  cache_a(2 / spec.r0 / spec.r0 * spec.v_rf ),
  cache_b(1.0/ spec.z0 / spec.z0 * spec.v_ec * spec.kappa),
  cache_c(2*M_PI*spec.omega) {
}

LQT::LQT() {}

Vector3D LQT::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  double a, b;
  double q = atom.GetCharge();

  a = cache_a * q * cos(cache_c*now);
  b = cache_b * q;

  Vector3D p(atom.GetPosition());
  return Vector3D((-a+b)*p[0], (a+b)*p[1], -2*b*p[2]);
}

// MathieuFunc::mathieu_param LQT::GetMathieuParam(const Util::ConstSIAtomProxy &atom) const {
//   double q, a;
//   q = 4 * atom.GetCharge() * spec.v_rf / (atom.GetMass() * spec.r0 * spec.r0 * spec.omega * spec.omega * 4 * M_PI * M_PI);
//   a = 4 * atom.GetCharge() * spec.v_ec / (atom.GetMass() * spec.z0 * spec.z0 * spec.omega * spec.omega * 4 * M_PI * M_PI);
  
//   return MathieuFunc::mathieu_param(q, a);
// }

// array<double, 3> LQT::GetSecularFrequency(const Util::ConstSIAtomProxy &atom) const {
//   array<double, 3> fsec;
  

//   return fsec;
// }
