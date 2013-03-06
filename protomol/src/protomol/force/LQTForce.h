#ifndef LQTFORCE_H
#define LQTFORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/output/LinearPaulTrap.h>
#include <protomol/output/LuaState.h>
#include <omp.h>

using std::string;
using LinearPaulTrap::Lqt;
using LuaState::LuaState;

namespace ProtoMol{
  template<class TBoundaryConditions>
  class LQTForce: public ExtendedForce {
  private:
    string filename;
    LuaState::LuaState L;
    Lqt trap;

    
  public:
    LQTForce(): 
      filename(""),
      L(),
      trap()
    {}

    LQTForce(const std::string& filename): 
      filename(filename), 
      L(filename),
      trap(L)
    {

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
    double force_conversion = Constant::SI::KCAL * Constant::SI::AVOGADRO * 1e-10;
    double position_conversion = 1e-10;
    double time_conversion = 1.0 / ProtoMol::Constant::SI::TIME_FS;

#pragma omp parallel for 
    for(unsigned int i=0;i<topo->atoms.size();i++)
      {
	Vector3D f;
	trap.GetForce((*positions)[i] * position_conversion, topo->time * time_conversion, f);
	(*forces)[i] += f * force_conversion;
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
    parameters.push_back(Parameter("-lqt_filename", Value(filename)));
  }
}
#endif
