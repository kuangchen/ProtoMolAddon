#ifndef _STRAYFIELD_H
#define _STRAYFIELD_H

#include <protomol/addon/Constants.h>
#include <protomol/addon/stray_field/StrayField.h>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <stdexcept>
#include <vector>
#include <string>

using std::string;
using namespace ProtoMolAddon::Constant;

using namespace ProtoMolAddon::StrayField;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

namespace ProtoMolAddon {
  namespace StrayField {
  
    template<class TBoundaryConditions>
      class StrayFieldForce: public ExtendedForce  {

    public:
      StrayFieldForce(): conf_fname(""), field() {}

      StrayFieldForce(const string& filename): 
	conf_fname(filename),
	field() {

	ifstream is(conf_fname);
	if (!is)
	  throw runtime_error(string("Cannot open file ") + conf_fname);

	is >> field;
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
      string conf_fname;
      StrayField field;
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
      for ( unsigned int i=0 ; i<topo->atoms.size() ; i++ )
	(*forces)[i] += field.GetForce(topo->atoms[i].scaledCharge * ToSI::charge) * ToSI::force;
	  //for ( unsigned int j=0 ; j<3 ; j++ )
	  //(*forces)[i][j] += field[j] * 1.60217646e-19 * conversion;
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
      parameters.push_back(Parameter("StrayFieldConf", Value(conf_fname)));
			 
    }

  }
}
#endif
