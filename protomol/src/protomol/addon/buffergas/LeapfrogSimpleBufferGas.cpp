#include <protomol/addon/buffergas/LeapfrogSimpleBufferGas.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/addon/util/SIAtomProxyArray.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/base/PMConstants.h>
#include <cmath>
#include <boost/algorithm/string.hpp>

namespace ProtoMolAddon {
  namespace BufferGas {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;
    
    LeapfrogSimpleBufferGas::NeutralAtom::NeutralAtom(const std::string &fname) 
      : rd(), engine(rd()) 
    {
      pt::ptree tree;
      pt::read_xml(fname, tree);
      
      m = tree.get<double>("ConfigRoot.BufferGas.m") * Constant::ToSI::mass;
      name = tree.get<std::string>("ConfigRoot.BufferGas.name");
      algorithm::trim(name);
      alpha = tree.get<double>("ConfigRoot.BufferGas.alpha");
      T = tree.get<double>("ConfigRoot.BufferGas.T");
      rho = tree.get<double>("ConfigRoot.BufferGas.rho");
      target_atom_name = tree.get<std::string>("ConfigRoot.BufferGas.target_atom_name");
      algorithm::trim(target_atom_name);
      
      sigma = sqrt(ProtoMol::Constant::SI::BOLTZMANN * T / m);
      vel_dist.param(std::normal_distribution<double>::param_type(0, sigma));
    }

    
    double LeapfrogSimpleBufferGas::NeutralAtom::GetCollisionInterval(Util::SIAtomProxy &ap) {
      double mu = m * ap.GetMass() / (m + ap.GetMass());
      double C4 = alpha/2 * 4.35974434e-18 * pow(5.2917721092e-11, 4);
      double gamma = 2 * M_PI * rho * sqrt(C4/mu);

      interval.param(std::exponential_distribution<double>::param_type(gamma));
      return interval(engine);
    }

    LeapfrogSimpleBufferGas::LeapfrogSimpleBufferGas() : neutral() {}

    LeapfrogSimpleBufferGas::LeapfrogSimpleBufferGas(const std::string &fname): 
      neutral(fname) {}

    void LeapfrogSimpleBufferGas::Initialize(ProtoMolApp *app) {
      ap_array_ptr.reset(new Util::SIAtomProxyArray(app));
    }

    void LeapfrogSimpleBufferGas::Update(double now, double dt) {
      for (auto &ap: *ap_array_ptr) 
	neutral.Collide(ap, dt);
    }

    void LeapfrogSimpleBufferGas::NeutralAtom::Collide(Util::SIAtomProxy &ap, double dt) {
      if (ap == target_atom_name || GetCollisionInterval(ap) > dt) return;

      Vector3D v(vel_dist(engine), vel_dist(engine), vel_dist(engine));
      const Vector3D vi = ap.GetVelocity();
      double mi = ap.GetMass();
      double v_rel_mag = (vi - v).norm();
      double m_total = m + mi;

      Vector3D v_com = (vi * mi + v * m) / m_total;
      Vector3D v_rel(uniform_dist(engine), uniform_dist(engine), uniform_dist(engine));
      v_rel.normalize();
      v_rel *= v_rel_mag;
      
      ap.SetVelocity(v_com + v_rel * m / m_total);
    }
  }
}
