#ifndef LQTWITHECFORCE_H
#define LQTWITHECFORCE_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/output/lqt_with_ec.h>

#include <omp.h>

namespace ProtoMol{
  using LQTWithEC::LQTWithEC;
  
  template<class TBoundaryConditions>
  class LQTWithECForce: public ExtendedForce
  {
  public:
    LQTWithECForce(): def_filename("")
    {}

    LQTWithECForce(const std::string& filename): def_filename(filename), lqt_with_ec(filename)
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
      return new LQTWithECForce(values[0]);
    };

    std::string def_filename;
    LQTWithEC lqt_with_ec;
  };
  
  template <class TBoundaryConditions>
  const std::string LQTWithECForce<TBoundaryConditions>::keyword("LQTWithEC");  

  template<class TBoundaryConditions>
  inline void LQTWithECForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double conversion = Constant::SI::KCAL * Constant::SI::AVOGADRO * 1e-10;

#pragma omp parallel for schedule(dynamic, 3)
    for( unsigned int i=0 ; i<topo->atoms.size() ; i++ ){
	Vector3D f;
	Vector3D pos((*positions)[i]);
	lqt_with_ec.GetForce((*positions)[i] * 1e-10, 
		     (*velocities)[i] * 1e-10 * ProtoMol::Constant::SI::TIME_FS/ ProtoMol::Constant::TIMEFACTOR , 
		     topo->time / ProtoMol::Constant::SI::TIME_FS, 
		     f);
	
	(*forces)[i] += f * conversion;
      }
  }


  template<class TBoundaryConditions>
  inline void LQTWithECForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void LQTWithECForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-lqt_with_ec_def_filename", Value(def_filename)));
  }
}
#endif
