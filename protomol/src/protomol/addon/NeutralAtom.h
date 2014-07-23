#ifndef __NEUTRAL_ATOM_H
#define __NEUTRAL_ATOM_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/LuaConfigReader.h>
#include <protomol/addon/ProtoMolIon.h>

#include <utility>
#include <iostream>
#include <random>

using namespace ProtoMol;
using namespace ProtoMolAddon::Lua;
using namespace std;


namespace ProtoMolAddon {
  class ProtoMolIon;

  class NeutralAtom {

  private:
    double mass;
    double polarizability;
    double temperature;
    double density;
    Vector3D position;
    Vector3D velocity;
    double sigma;

    typedef std::piecewise_constant_distribution<double> pc_dist;
    typedef pc_dist::param_type pc_dist_param_type;

    vector<double> vn_edge;
    vector<double> vn_center;
    vector<double> vn_weight;

    vector<double> theta_edge;
    vector<double> theta_center;
    vector<double> theta_weight;

    std::random_device rd;
    std::uniform_real_distribution<double> phi_dice;
    std::uniform_real_distribution<double> collision_dice;

    std::normal_distribution<double> vn_dice;
    

    typedef map<double, double> scalar_lktb;
    scalar_lktb f_lktb;

    typedef map<double, pc_dist> pc_dist_lktb;
    pc_dist_lktb diff_cross_section_lktb;

    
    bool ShouldCollide(const ProtoMolIon &ion, double dt);
    bool ShouldCollide2(const ProtoMolIon &ion, double dt);
    void Collide(ProtoMolIon &ion);
    void SampleVelocity(const ProtoMolIon &ion);
    void SampleVelocity2(const ProtoMolIon &ion);
    double GetThermalAverageCollisionRate(const ProtoMolIon &ion);
    double GetCollisionRate(const ProtoMolIon &ion);

    void LoadFLktb();
    void LoadThetaCache(unsigned int count=100);
    void LoadVnCache(unsigned int count=100);
    void LoadCrossSection(LuaConfigReader &reader);
    pc_dist &SelectThetaDist(double energy_in_kelvin);

  public:
    NeutralAtom(double m=0, double pl=0, double t=0, double d=0, const Vector3D &pos=Vector3D(), const Vector3D &vel=Vector3D());

    NeutralAtom(LuaConfigReader &reader);
    void CollideAll(ProtoMolApp *app, double dt);
    void CollideAll2(ProtoMolApp *app, double dt);
    void SampleVelocityTest(const string &fname);
    void ScatteringTest(const string &fname);
    

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
