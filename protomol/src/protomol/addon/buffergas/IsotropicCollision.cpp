#include <protomol/addon/buffergas/IsotropicCollision.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <protomol/base/PMConstants.h>
#include <cmath>
#include <boost/algorithm/string.hpp>

namespace pt = boost::property_tree;
namespace algorithm = boost::algorithm;

namespace ProtoMolAddon {
  namespace BufferGas {

    IsotropicCollision::NeutralAtom::NeutralAtom(const std::string &fname) 
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
      target = tree.get<std::string>("ConfigRoot.BufferGas.target");
      algorithm::trim(target);
      
      sigma = sqrt(ProtoMol::Constant::SI::BOLTZMANN * T / m);
      vel_dist.param(std::normal_distribution<double>::param_type(0, sigma));
    }
    
    double IsotropicCollision::NeutralAtom::GetTimeToNextCollision(Util::SIAtomProxy &ap) {
      double mu = m * ap.GetMass() / (m + ap.GetMass());
      double C4 = alpha/2 * 4.35974434e-18 * pow(5.2917721092e-11, 4);
      double gamma = 2 * M_PI * rho * sqrt(C4/mu);
      //std::cout << "gamma = " << gamma << std::endl;
      time_to_next_collision.param(std::exponential_distribution<double>::param_type(gamma));

      return time_to_next_collision(engine);
    }

    IsotropicCollision::IsotropicCollision() : Collision(), neutral() {}

    IsotropicCollision::IsotropicCollision(const std::string &fname): 
      Collision(), 
      neutral(fname) {}

    void IsotropicCollision::CollideEach(Util::SIAtomProxy &ap, double dt) {
      neutral.Collide(ap, dt);
    }

    void IsotropicCollision::NeutralAtom::Collide(Util::SIAtomProxy &ap, double dt) {
      if (ap.GetName() != target || GetTimeToNextCollision(ap) > dt) 
	return;

      v = ProtoMol::Vector3D(vel_dist(engine), vel_dist(engine), vel_dist(engine));

      double v_rel_mag = (ap.GetVelocity() - v).norm();
      double m_total = m + ap.GetMass();

      ProtoMol::Vector3D v_com = (ap.GetVelocity() * ap.GetMass() + v * m) / m_total;
      ProtoMol::Vector3D v_rel(uniform_dist(engine), uniform_dist(engine), uniform_dist(engine));
      v_rel.normalize();
      v_rel *= v_rel_mag;
      
      ap.SetVelocity(v_com + v_rel * m / m_total);
    }
  }
}
