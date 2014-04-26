#ifndef _SIMPLE_DAMPING_FORCE_H
#define _SIMPLE_DAMPING_FORCE_H

#include <string>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/SIAtomProxy.h>

namespace ProtoMolAddon {
  namespace Damping {

    using namespace ProtoMolAddon::Constant;

    template<class TBoundaryConditions, class Damping>
    class GenericDampingForce: public ExtendedForce {
    private: 
      std::string fname;
      Damping damping;

    public:
      GenericDampingForce() {}
      GenericDampingForce(const string& fname):
	fname(fname),
	damping(fname)
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
      virtual std::string getKeyword() const { return keyword; }
      virtual std::string getIdNoAlias() const { return keyword; }

    private:
      virtual Force* doMake(const std::vector<Value> &values) const
	{
	  return new GenericDampingForce(values[0]);
	};

    };
  
    template <class TBoundaryConditions, class Damping>
    const std::string GenericDampingForce<TBoundaryConditions, Damping>::keyword(Damping::GetName());  

    template<class TBoundaryConditions, class Damping>
    inline void GenericDampingForce<TBoundaryConditions, Damping>::evaluate(const GenericTopology* topo,
									   const Vector3DBlock* positions,
									   const Vector3DBlock *velocities,
									   Vector3DBlock* forces,
									   ScalarStructure* energies)
    {
      for(unsigned int i=0;i<topo->atoms.size();i++) {
	Util::ConstSIAtomProxy atom(topo, positions, velocities, i);
	(*forces)[i] += damping.GetForce(atom, topo->time * ToSI::time) / ToSI::force;
      }
    }


    template<class TBoundaryConditions, class Damping>
    inline void GenericDampingForce<TBoundaryConditions, Damping>::parallelEvaluate(const GenericTopology* topo,
										   const Vector3DBlock* positions,
										   const Vector3DBlock *velocities,
										   Vector3DBlock* forces,
										   ScalarStructure* energies)
    {
      evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
    }

    template<class TBoundaryConditions, class Damping>
    inline void GenericDampingForce<TBoundaryConditions, Damping>::getParameters(std::vector<Parameter>& parameters) const
    {
      parameters.push_back(Parameter("-simple_damping_spec", Value(fname)));
    }
  }
}

#endif

