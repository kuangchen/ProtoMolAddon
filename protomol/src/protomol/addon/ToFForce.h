#ifndef _TOF_FORCE_H
#define _TOF_FORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/addon/ToF.h>
#include <protomol/addon/Constants.h>
#include <omp.h>

using namespace std;
using namespace ProtoMol::Constant;
using namespace ProtoMol::Constant::SI;
using namespace ProtoMolAddon::Constant;

namespace ProtoMolAddon{

  template<class TBoundaryConditions>
  class ToFForce: public ExtendedForce {
  private: 
    string def;
    ToF tof;

  public:
    ToFForce() {}
    ToFForce(const string& def):
      def(def),
      tof(def)
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
      return new ToFForce(values[0]);
    };

  };
  
  template <class TBoundaryConditions>
  const std::string ToFForce<TBoundaryConditions>::keyword("ToFForce");  

  template<class TBoundaryConditions>
  inline void ToFForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    for(unsigned int i=0;i<topo->atoms.size();i++) {
	Vector3D f;
	Vector3D pos((*positions)[i]);
	tof.GetForce(topo->atoms[i].scaledCharge * CHARGE_CONV * ELECTRON_CHARGE, 
		     pos * POSITION_CONV, 
		     topo->time * TIME_CONV, 
		     f);
	//	cout << scientific << "i = " << i << "\tf = " << f << "\n";
	(*forces)[i] += f * FORCE_CONV;
      }
  }


  template<class TBoundaryConditions>
  inline void ToFForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void ToFForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-tof_def", Value(def)));
			 
  }
}

#endif
