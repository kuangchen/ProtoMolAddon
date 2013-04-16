#ifndef __NEUTRAL_ATOM_H
#define __NEUTRAL_ATOM_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/LuaConfigReader.h>

#include <utility>
#include <iostream>
#include <random>

using namespace ProtoMol;
using namespace ProtoMolAddon::Lua;
using namespace std;

namespace ProtoMolAddon {

  class NeutralAtom {
  private:
    double mass;
    double polarizability;
    double temperature;
    double density;
    Vector3D position;
    Vector3D velocity;
    
    double sigma;
    vector<double> r_cache;
    vector<double> r_cache_middle;
    vector<double> w_r_cache;
    vector<double> theta_cache;
    vector<double> theta_cache_middle;
    vector<double> w_theta_cache;

    std::random_device rd;
    std::uniform_real_distribution<double> uniform_dist;
    std::piecewise_constant_distribution<double> mag_dist;
    std::piecewise_constant_distribution<double> polar_angle_dist;
    std::piecewise_constant_distribution<double> polar_angle_dist2;
    std::uniform_real_distribution<double> azimuthal_angle_dist;


    std::map<double, double> lookup_table1;

    void Collide(double m, Vector3D &v);
    Vector3D SampleVelocity(const Vector3D &v);
    double GetAverageCollisionRate(double m, double q, const Vector3D &v);
    double GetCollisionRate(double m, double q, const Vector3D &v);
    double GetLangevinRate(double m, double q, const Vector3D &v);

  public:
    NeutralAtom(double m=0, double pl=0, double t=0, double d=0, const Vector3D &pos=Vector3D(), const Vector3D &vel=Vector3D());

    NeutralAtom(LuaConfigReader &reader);
    void CollideAll(ProtoMolApp *app, double dt);
    void LoadLookupTable();
    void LoadCache(); 
    void Test();

    friend std::ostream& operator<< (std::ostream& os, const NeutralAtom &neutral){
      os << "Dumping Neutral Information ... " 
	 << "\n---------------------------------- " 
	 << "\nNeutral: Mass           = " << neutral.mass 
	 << "\n       : Polarizability = " << neutral.polarizability 
	 << "\n       : Temperature    = " << neutral.temperature 
	 << "\n       : Density        = " << neutral.density
	 << "\n       : C4 / Q^2       = " << neutral.polarizability * 1.6e-19 * 1.6e-19 / (4*M_PI*8.85e-12) << "\n";
  
      return os;
    }
  

  };
}

#endif
