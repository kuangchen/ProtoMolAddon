#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/type/Vector3D.h>

#include <stdexcept>
#include <iostream>
#include <boost/program_options.hpp>


using std::istream;
using std::ostream;
using namespace ProtoMolAddon::IonTrap;
namespace po = boost::program_options;


HarmonicTrap::HarmonicTrapSpec::HarmonicTrapSpec(const std::string &fname) {
  std::ifstream is(fname.c_str());

  if (!is)
    throw std::runtime_error(string("Cannot open file ") + fname);

  po::variables_map vm; 
  po::options_description desc("HarmonicTrap Spec"); 
  desc.add_options()
    ("HarmonicTrap.omega_x", po::value<double>(&omega[0])->required(), "trap frequency x")
    ("HarmonicTrap.omega_y", po::value<double>(&omega[1])->required(), "trap frequency y")
    ("HarmonicTrap.omega_z", po::value<double>(&omega[2])->required(), "trap frequency z");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);
}

HarmonicTrap::HarmonicTrap() {}

HarmonicTrap::HarmonicTrap(const HarmonicTrap::HarmonicTrapSpec &spec) :
  spec(spec) {
}


Vector3D HarmonicTrap::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  Vector3D f;
  double m = atom.GetMass();
  Vector3D p(atom.GetPosition());

  for (int i=0; i<3; i++)
    f[i] = spec.omega[i] * spec.omega[i] * p[i];

  return f * (-m);
}


