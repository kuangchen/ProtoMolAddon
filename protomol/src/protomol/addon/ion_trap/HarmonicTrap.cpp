#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <protomol/addon/util/ConstSIAtomProxy.h>
#include <protomol/type/Vector3D.h>

#include <stdexcept>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>


namespace ProtoMolAddon {
  namespace IonTrap {

    namespace pt = boost::property_tree;
    
    HarmonicTrap::Spec::Spec(const std::string &fname) {
      pt::ptree tree;
      pt::read_xml(fname, tree);

      omega = tree.get<Vector3D>("ConfigRoot.HarmonicTrap.omega");
    }

    HarmonicTrap::HarmonicTrap(const HarmonicTrap::Spec &spec) :
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


  }
}
