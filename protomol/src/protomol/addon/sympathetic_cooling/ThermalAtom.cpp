#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <cmath>
#include <string>
#include <boost/program_options.hpp>
#include <fstream>

namespace po = boost::program_options;

using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::SympatheticCooling;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

ThermalAtom::ThermalAtom() {}

ThermalAtom::ThermalAtom(const string &fname) :
  position(0, 0, 0),
  rd(), 
  dice(0, sqrt(SI::BOLTZMANN*temperature/mass))
{
  std::ifstream is(fname.c_str());

  if (!is)
    throw runtime_error(string("Cannot open file " + fname));

  po::variables_map vm; 
  po::options_description desc("ThermalAtom spec"); 
  desc.add_options()
    ("ThermalAtom.mass", po::value<double>(&mass)->required(), "atom mass [amu]")
    ("ThermalAtom.density", po::value<double>(&density)->required(), "atom density [m^(-3)]")
    ("ThermalAtom.temperature", po::value<double>(&temperature)->required(), "atom temperature [K]")
    ("ThermalAtom.polarizability", po::value<double>(&polarizability)->required(), "atomic polarizability [cm^3]");

  if (mass<0 || density<0 || temperature<0) 
    throw runtime_error("invalid input for mass/density/temperature");

  po::store(po::parse_config_file(is, desc, true), vm); 
  po::notify(vm);
  mass = mass * ToSI::mass;
  
  C4 = polarizability * ToSI::charge * ToSI::charge / (4*M_PI*epsilon_0);

}

void ThermalAtom::Resample() {
  velocity = Vector3D(dice(rd), dice(rd), dice(rd));
}

