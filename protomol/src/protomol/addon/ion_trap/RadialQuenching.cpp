#include <protomol/addon/ion_trap/RadialQuenching.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/addon/util/SIAtomProxy.h>

using namespace ProtoMolAddon::IonTrap;

namespace pt = boost::property_tree;

RadialQuenching::Spec::Spec(const std::string &fname) {
  pt::ptree tree;
  pt::read_xml(fname, tree);

  label = tree.get<std::string>("RadialQuenchingSpec.label");
  r0 = tree.get<double>("RadialQuenchingSpec.r0");
  v = tree.get<double>("RadialQuenchingSpec.v");
}

RadialQuenching::RadialQuenching(const RadialQuenching::Spec &spec) :
  spec(spec) {}

Vector3D RadialQuenching::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  Vector3D p(atom.GetPosition());

  return Vector3D(-2*p[0], 2*p[1], 0) * (spec.v / spec.r0 / spec.r0 * atom.GetCharge());
}
