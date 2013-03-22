#ifndef _DAMPING_FORCE_H
#define _DAMPING_FORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/addon/Damping.h>

using namespace std;
using namespace ProtoMol::Constant;

namespace ProtoMolAddon{

  template<class TBoundaryConditions>
  class DampingForce: public ExtendedForce {
  private: 
    string def;
    Damping d;

  public:
    DampingForce() {}
    DampingForce(const string& def):
      def(def),
      d(def)
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
      return new DampingForce(values[0]);
    };

  };
  
  template <class TBoundaryConditions>
  const std::string DampingForce<TBoundaryConditions>::keyword("DampingForce");  

  template<class TBoundaryConditions>
  inline void DampingForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double force_conversion = SI::KCAL * SI::AVOGADRO * 1e-10;
    double position_conversion = 1e-10;
    double time = topo->time/SI::TIME_FS;
    double velocity_conversion = 1e-10 * SI::TIME_FS / TIMEFACTOR;
    for(unsigned int i=0;i<topo->atoms.size();i++)
      {
	Vector3D f;
	Vector3D vel((*velocities)[i]); 
	d.GetForce(vel * velocity_conversion, time, f);
	(*forces)[i] += f * force_conversion;
      }
  }


  template<class TBoundaryConditions>
  inline void DampingForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void DampingForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-damping_def", Value(def)));
  }
}

#endif
