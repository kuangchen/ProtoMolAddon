#ifndef _SEGMENTED_TRAP_FORCE_H
#define _SEGMENTED_TRAP_FORCE_H

#include <protomol/addon/segmented_trap/NewtonDiff.h>
#include <protomol/addon/segmented_trap/DummyElectrode.h>
#include <protomol/addon/segmented_trap/SegmentedTrap.h>
#include <protomol/addon/Constants.h>
#include <protomol/type/Vector3D.h>
#include <vector>

using namespace ProtoMolAddon::Constant;
using namespace ProtoMol;

namespace ProtoMolAddon {
  namespace SegmentedTrap {

    using namespace std;

    template<class TBoundaryConditions>
    class SegmentedTrapForce: public ExtendedForce  {
    private:
      string conf_fname;
      SegmentedTrap<DummyElectrode, NewtonDiff<DummyElectrode> > field;

    public:
      SegmentedTrapForce(): conf_fname(""), field() {}

      SegmentedTrapForce(const string& filename): 
	conf_fname(filename),
	field() {

	ifstream is(conf_fname);
	if (!is)
	  throw runtime_error(string("Cannot open file ") + conf_fname);

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
	  return new SegmentedTrapForce(values[0]);
	};

    };
  
    template <class TBoundaryConditions>
    const std::string SegmentedTrapForce<TBoundaryConditions>::keyword("SegmentedTrap");  

    template<class TBoundaryConditions>
    inline void SegmentedTrapForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
								  const Vector3DBlock* positions,
								  const Vector3DBlock *velocities,
								  Vector3DBlock* forces,
								  ScalarStructure* energies)
    {
      Vector3D f;
      Vector3D x;
      for (int i=0;i<topo->atoms.size();i++) {
	x = (*positions)[i];
	field.GetForce(topo->atoms[i].scaledCharge * ToSI::charge, x * ToSI::position,0 , f);
	(*forces)[i] += f;
      }
    }


    template<class TBoundaryConditions>
    inline void SegmentedTrapForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
									  const Vector3DBlock* positions,
									  const Vector3DBlock *velocities,
									  Vector3DBlock* forces,
									  ScalarStructure* energies)
    {
      evaluate(topo, positions, velocities, forces, energies); // Not implemented right now
    }

    template<class TBoundaryConditions>
    inline void SegmentedTrapForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
    {
      parameters.push_back(Parameter("-SegmentedTrapConf", Value(conf_fname)));
			 
    }

  }
}

#endif
