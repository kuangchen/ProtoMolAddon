#ifndef _HARMONIC_TRAP_FORCE_H_
#define _HARMONIC_TRAP_FORCE_H_

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <fstream>
#include <string>

using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace IonTrap {
    
    using namespace ProtoMolAddon::Util;
    
    template<class TBoundaryConditions>
    class HarmonicTrapForce: public ExtendedForce {
    private: 
      std::string fname;
      HarmonicTrap trap;

    public:
      HarmonicTrapForce(): fname(""), trap() {}

      HarmonicTrapForce(const string& fname): fname(fname), trap(fname) {}
    
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
      for(unsigned int i=0;i<topo->atoms.size();i++) {
	ConstSIAtomProxy atom(topo, positions, velocities, i);
	(*forces)[i] += trap.GetForce(atom, topo->time * ToSI::time) / ToSI::force;
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
      parameters.push_back(Parameter("-ht_spec", Value(fname)));
			 
    }

  }
}

#endif
