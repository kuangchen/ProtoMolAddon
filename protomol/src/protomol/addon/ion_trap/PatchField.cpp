#include <protomol/addon/ion_trap/PatchField.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/addon/util/ConstSIAtomProxy.h>

namespace ProtoMolAddon {
  namespace IonTrap {

    namespace pt = boost::property_tree;

    PatchField::Spec::Spec(const std::string &fname) {
      pt::ptree tree;
      pt::read_xml(fname, tree);

      label = tree.get<std::string>("ConfigRoot.PatchFieldSpec.label");
      r0 = tree.get<double>("ConfigRoot.PatchFieldSpec.r0");
      v = tree.get<double>("ConfigRoot.PatchFieldSpec.v");
    }

    PatchField::PatchField(const PatchField::Spec &spec) : spec(spec) {}

    Vector3D PatchField::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
      const Vector3D &p = atom.GetPosition();
      double c(atom.GetCharge());
      return Vector3D(-2*p[0], 2*p[1], 0) * (spec.v / spec.r0 / spec.r0 * c);
    }


  }
}
