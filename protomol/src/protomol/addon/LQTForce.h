#ifndef LQTFORCE_H
#define LQTFORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>

#include <protomol/addon/LinearPaulTrap.h>
#include <protomol/addon/LuaConfigReader.h>

#include <vector>
#include <string>

using namespace ProtoMolAddon::LinearPaulTrap;
using namespace ProtoMolAddon::Lua;
using namespace ProtoMol::Constant;


namespace ProtoMolAddon {

  template<class TBoundaryConditions>
  class LQTForce: public ExtendedForce {
  private:
    string filename;
    LuaConfigReader reader;
    Lqt trap;
    
    static const double force_conv;
    static const double position_conv;
    static const double time_conv;
    static const double charge_conv;

  public:
    LQTForce(): 
      filename(""),
      reader(),
      trap()
    {}

    LQTForce(const std::string& filename): 
      filename(filename), 
      reader(filename),
      trap(reader)
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

    for(unsigned int i=0;i<topo->atoms.size();i++) 
      (*forces)[i] += trap.GetForce(topo->atoms[i].scaledCharge * charge_conv, 
				    (*positions)[i] * position_conv, 
				    topo->time * time_conv) * force_conv;
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

  template<class TBoundaryConditions> const double LQTForce<TBoundaryConditions>::force_conv = SI::KCAL * SI::AVOGADRO * 1e-10;
  template<class TBoundaryConditions> const double LQTForce<TBoundaryConditions>::position_conv = 1e-10;
  template<class TBoundaryConditions> const double LQTForce<TBoundaryConditions>::time_conv = 1.0 / SI::TIME_FS;
  template<class TBoundaryConditions> const double LQTForce<TBoundaryConditions>::charge_conv = 1.0 / SQRTCOULOMBCONSTANT;

  
}



#endif
