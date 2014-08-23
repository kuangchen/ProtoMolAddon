#ifndef __ISOTROPIC_COLLISION_H_
#define __ISOTROPIC_COLLISION_H_

#include <protomol/addon/buffergas/Collision.h>
#include <protomol/type/Vector3D.h>
#include <random>
#include <string>

namespace ProtoMolAddon {
  namespace BufferGas {

    class IsotropicCollision : public Collision {

      struct NeutralAtom {
      private:
	double m;
	std::string name;
	double alpha;
	double T;
	double rho;
	std::string target;

	std::random_device rd;
	std::default_random_engine engine;
	std::normal_distribution<double> vel_dist;
	std::exponential_distribution<double> time_to_next_collision;
	std::uniform_real_distribution<double> uniform_dist;
	double sigma;
	ProtoMol::Vector3D v;

      private:
	double GetTimeToNextCollision(Util::SIAtomProxy &ap);

      public:
	NeutralAtom() {};
	NeutralAtom(const std::string &fname);
	void Collide(Util::SIAtomProxy &ap, double dt);
      };

    private:
      NeutralAtom neutral;

    protected:
      void CollideEach(Util::SIAtomProxy &ap, double dt);

    public:
      static const std::string GetKeyword() { return "Isotropic"; };
      IsotropicCollision();
      IsotropicCollision(const std::string &fname);
    };

  }
}


#endif
