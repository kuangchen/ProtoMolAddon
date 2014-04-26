#ifndef _LQT_FORCE_H
#define _LQT_FORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/addon/ion_trap/LQT.h>
#include <string>

using namespace ProtoMol::Constant;

namespace ProtoMolAddon {
  namespace IonTrap {

    //using namespace ProtoMolAddon::Util;

    template<class TBoundaryConditions>
    class LQTForce: public ExtendedForce {
    private:
      std::string spec_fname;
      LQT trap;

    public:
      LQTForce() {}

      LQTForce(const std::string& fname): 
	spec_fname(fname), trap(fname)
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
      virtual Force* doMake(const std::vector<Value> &values) const {
	return new LQTForce(values[0]);
      };

    };
  
    template <class TBoundaryConditions>
    const std::string LQTForce<TBoundaryConditions>::keyword("LQT");  

    template<class TBoundaryConditions>
    inline void LQTForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
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
    inline void LQTForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
								const Vector3DBlock* positions,
								const Vector3DBlock *velocities,
								Vector3DBlock* forces,
								ScalarStructure* energies)
    {
      evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
    }

    template<class TBoundaryConditions>
    inline void LQTForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
    {
      parameters.push_back(Parameter("-lqt_spec", Value(spec_fname)));
    }

  }
}



#endif
