#include <fstream>
#include <cmath>
#include <boost/program_options.hpp>
#include <protomol/addon/ion_trap/LQT.h>
#include <protomol/addon/util/SIAtomProxy.h>

namespace po = boost::program_options;
using namespace ProtoMolAddon::IonTrap;

LQT::LQTSpec::LQTSpec(const string &fname) {
  std::ifstream is(fname.c_str());

  if (!is)
    throw runtime_error(string("Cannot open file ") + fname);

  po::variables_map vm; 
  po::options_description desc("LQT Spec"); 
  desc.add_options()
    ("LQT.r0", po::value<double>(&r0)->required(), "trap radius")
    ("LQT.z0", po::value<double>(&z0)->required(), "trap length")
    ("LQT.omega", po::value<double>(&omega)->required(), "trap frequency")
    ("LQT.v_rf", po::value<double>(&v_rf)->required(), "rf voltage (amplitude)")
    ("LQT.v_ec", po::value<double>(&v_ec)->required(), "end-cap voltage")
    ("LQT.kappa", po::value<double>(&kappa)->required(), "geometrical factor");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);
}


LQT::LQT(const LQT::LQTSpec &spec) :
  spec(spec), 
  cache_a(2 / spec.r0 / spec.r0 * spec.v_rf ),
  cache_b(1.0/ spec.z0 / spec.z0 * spec.v_ec * spec.kappa),
  cache_c(2*M_PI*spec.omega)
{
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
