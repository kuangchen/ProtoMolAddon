#ifndef _HARMONIC_TRAP_FORCE_H_
#define _HARMONIC_TRAP_FORCE_H_

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <protomol/addon/Constants.h>
#include <fstream>
#include <string>

using std::ifstream;
using std::string;
using namespace ProtoMol;
using namespace std;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon;
using ProtoMolAddon::IonTrap::HarmonicTrap;

namespace ProtoMolAddon {
  namespace IonTrap {
  
    template<class TBoundaryConditions>
    class HarmonicTrapForce: public ExtendedForce {
    private: 
      string conf;
      HarmonicTrap trap;

    public:
      HarmonicTrapForce(): conf(""), trap() {}
      HarmonicTrapForce(const string& conf): conf(conf) {
	ifstream is(conf);
	is >> trap;
      }
    
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
      for(unsigned int i=0;i<topo->atoms.size();i++) 
	(*forces)[i] += trap.GetForce((*positions)[i] * Constant::position_conv, 
				      topo->atoms[i].scaledMass * SI::AMU) * Constant::force_conv;
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
      parameters.push_back(Parameter("-conf", Value(conf)));
			 
    }

  }
}

#endif
