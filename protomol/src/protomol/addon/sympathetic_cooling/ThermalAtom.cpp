#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <stdexcept>
#include <cmath>
#include <boost/program_options.hpp>
#include <fstream>

namespace po = boost::program_options;

using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::SympatheticCooling;
using namespace ProtoMol::Constant;

ThermalAtom::ThermalAtomSpec::ThermalAtomSpec() {}

ThermalAtom::ThermalAtomSpec::ThermalAtomSpec(const std::string &fname) {
  std::ifstream is(fname.c_str());
  
  if (!is)
    throw std::runtime_error(string("Cannot open file " + fname));

  po::variables_map vm; 
  po::options_description desc("ThermalAtom spec"); 
  desc.add_options()
    ("ThermalAtom.name", po::value<std::string>(&name)->required(), "atom name")
    ("ThermalAtom.mass", po::value<double>(&mass)->required(), "atom mass [amu]")
    ("ThermalAtom.density", po::value<double>(&density)->required(), "atom density [m^(-3)]")
    ("ThermalAtom.temperature", po::value<double>(&temperature)->required(), "atom temperature [K]")
    ("ThermalAtom.polarizability", po::value<double>(&polarizability)->required(), "atomic polarizability [cm^3]");

  if (mass<0 || density<0 || temperature<0) 
    throw std::runtime_error("invalid input for mass/density/temperature");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);

  mass = mass * ToSI::mass;
}


ThermalAtom::ThermalAtom() {}

ThermalAtom::ThermalAtom(const ThermalAtom::ThermalAtomSpec &spec) :
  spec(spec),
  C4(spec.polarizability * ToSI::charge * ToSI::charge / (4*M_PI*epsilon_0)),
  position(0, 0, 0),
  velocity(0, 0, 0),
  rd(), 
  dice(0, sqrt(SI::BOLTZMANN*spec.temperature/spec.mass))
{
}

void ThermalAtom::Resample() {
  velocity = Vector3D(dice(rd), dice(rd), dice(rd));
}

