#ifndef STRAYFIELD_H
#define STRAYFIELD_H

#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>

#include <protomol/addon/LuaConfigReader.h>
#include <vector>
#include <string>

using namespace std;
using namespace ProtoMolAddon::Lua;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

namespace ProtoMolAddon {
  
  template<class TBoundaryConditions>
  class StrayFieldForce: public ExtendedForce  {

  public:
    StrayFieldForce(): 
    sf_def_filename(""), 
      reader(),
      field(0)
    {}

    StrayFieldForce(const string& filename): 
      sf_def_filename(filename),
      reader(filename),
      field(reader.GetValue<vector<double> >("field"))
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
      return new StrayFieldForce(values[0]);
    };

    string sf_def_filename;
    LuaConfigReader reader;
    vector<double> field;
  };
  
  template <class TBoundaryConditions>
  const std::string StrayFieldForce<TBoundaryConditions>::keyword("StrayFieldForce");  

  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions,
						      const Vector3DBlock *velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* energies)
  {
    double conversion = SI::KCAL * SI::AVOGADRO * 1e-10;

    for ( unsigned int i=0 ; i<topo->atoms.size() ; i++ )
      for ( unsigned int j=0 ; j<3 ; j++ )
	(*forces)[i][j] += field[j] * 1.60217646e-19 * conversion;
  }


  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
							      const Vector3DBlock* positions,
							      const Vector3DBlock *velocities,
							      Vector3DBlock* forces,
							      ScalarStructure* energies)
  {
    evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
  }

  template<class TBoundaryConditions>
  inline void StrayFieldForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
  {
    parameters.push_back(Parameter("-sf_def_filename", Value(sf_def_filename)));
			 
  }
}
#endif
