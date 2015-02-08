#ifndef _FORCE_TEMPLATE_H
#define _FORCE_TEMPLATE_H

#include <string>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/ConstSIAtomProxy.h>

namespace ProtoMolAddon {
  namespace Template {

    using namespace ProtoMol;
    using namespace ProtoMolAddon::Constant;

    template<class TBoundaryConditions, class F>
    class GenericForce: public ExtendedForce {
    private: 
      std::string fname;
      F f;

    public:
      GenericForce() {}
      GenericForce(const std::string& fname):
	fname(fname),
	f(fname)
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
	  return new GenericForce(values[0]);
	};

    };
  
    template <class TBoundaryConditions, class F>
    const std::string GenericForce<TBoundaryConditions, F>::keyword(F::GetName());  

    template<class TBoundaryConditions, class F>
    inline void GenericForce<TBoundaryConditions, F>::evaluate(const GenericTopology* topo,
									   const Vector3DBlock* positions,
									   const Vector3DBlock *velocities,
									   Vector3DBlock* forces,
									   ScalarStructure* energies)
    {
      for(unsigned int i=0;i<topo->atoms.size();i++) {
	Util::ConstSIAtomProxy ap(topo, positions, velocities, i);
	(*forces)[i] += f.GetForce(ap, topo->time * ToSI::time) / ToSI::force;
      }
    }


    template<class TBoundaryConditions, class F>
    inline void GenericForce<TBoundaryConditions, F>::parallelEvaluate(const GenericTopology* topo,
										   const Vector3DBlock* positions,
										   const Vector3DBlock *velocities,
										   Vector3DBlock* forces,
										   ScalarStructure* energies)
    {
      evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
    }

    template<class TBoundaryConditions, class F>
    inline void GenericForce<TBoundaryConditions, F>::getParameters(std::vector<Parameter>& parameters) const
    {
      parameters.push_back(Parameter(F::GetParameterName(), Value(fname)));
    }
  }
}

#endif 
