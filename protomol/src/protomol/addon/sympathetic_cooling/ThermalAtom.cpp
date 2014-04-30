#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace pt = boost::property_tree;

using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::SympatheticCooling;
using namespace ProtoMol::Constant;

ThermalAtom::Spec::Spec() {}

ThermalAtom::Spec::Spec(const std::string &fname) {
  pt::ptree tree;
  pt::read_xml(fname, tree);

  name = tree.get<std::string>("Atom.name");
  mass = tree.get<double>("Atom.mass") * ToSI::mass;
  density = tree.get<double>("Atom.density");
  temperature = tree.get<double>("Atom.temperature");
  polarizability = tree.get<double>("Atom.polarizability");

  if (mass<0 || density<0 || temperature<0) 
    throw std::runtime_error("invalid input for mass/density/temperature");
}

ThermalAtom::ThermalAtom() {}

ThermalAtom::ThermalAtom(const ThermalAtom::Spec &spec) :
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

