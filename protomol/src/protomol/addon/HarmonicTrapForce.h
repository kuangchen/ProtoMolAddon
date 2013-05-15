#ifndef _HARMONIC_TRAP_FORCE_H
#define _HARMONIC_TRAP_FORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/addon/HarmonicTrap.h>
#include <protomol/addon/Constants.h>

using namespace std;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Constant;

namespace ProtoMolAddon{

  template<class TBoundaryConditions>
  class HarmonicTrapForce: public ExtendedForce {
  private: 
    string def;
    HarmonicTrap trap;

  public:
    HarmonicTrapForce() {}
    HarmonicTrapForce(const string& def):
      def(def),
      trap(def)
    {}
    
    virtual void evaluate(const GenericTopology* topo,
			  const Vector3DBlock* positions,
                          const Vector3DBlock *velocities,
			  Vector3DBlock* forces,
			  ScalarStructure* energies);

    virtual void parallelEvaluate(const GenericTopology* topo,
				  const Vector3DBlock* positions,
				  const Vector3DBlock *velocities,
				  Vector3DBlock* forces,
				  ScalarStructure* energies);

    virtual void getParameters(std::vector<Parameter>& ) const; 
    static const std::string keyword;
    virtual std::string getKeyword() const{return keyword;}
    virtual std::string getIdNoAlias() const{return keyword;}

  private:
    virtual Force* doMake(const std::vector<Value> &values) const
    {
      return new HarmonicTrapForce(values[0]);
    };

  };
  
  template <class TBoundaryConditions>
  const std::string HarmonicTrapForce<TBoundaryConditions>::keyword("HarmonicTrapForce");  

  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double force_conversion = SI::KCAL * SI::AVOGADRO * 1e-10;

    for(unsigned int i=0;i<topo->atoms.size();i++) {
	Vector3D f;
	Vector3D pos((*positions)[i]);
	double mass = topo->atoms[0].scaledMass * SI::AMU;
	trap.GetForce(pos * POSITION_CONV, mass, f);
	(*forces)[i] += f * force_conversion;
      }
  }


  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void HarmonicTrapForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-ht_def", Value(def)));
			 
  }
}

#endif
