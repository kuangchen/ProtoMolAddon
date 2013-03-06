%module MultiGridForce
%{
#include <protomol/base/Report.h>
#include <protomol/config/Value.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomPairTwo.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/switch/UniversalSwitchingFunction.h>
#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/ComplementSwitchingFunction.h>
#include <protomol/switch/CutoffSwitchingFunction.h>
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
//#include <protomol/force/coulomb/CoulombBornRadiiForce.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/table/LennardJonesTableForce.h>
#include <protomol/force/table/CoulombTableForce.h>
#include <protomol/topology/CellListEnumeratorStandard.h>
#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include "Hermite.h"
#include "NonbondedMultiGridSystemForce.h"

#include <protomol/base/Report.h>
using namespace ProtoMol;
%}

%include "std_string.i"
%include <protomol/type/Real.h>
//%include <protomol/config/Value.h>
%include <protomol/config/Parameter.h>
%include <protomol/type/Vector3D.h>
%include <interface/type/Vector3DBlock.i>
%include <protomol/topology/Angle.h>
%include <protomol/topology/Atom.h>
%include <protomol/topology/AtomType.h>
%include <protomol/topology/Bond.h>
%include <protomol/topology/Torsion.h>
%include <protomol/topology/Molecule.h>
%include <protomol/type/ScalarStructure.h>
%include <protomol/topology/GenericTopology.h>
%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/base/Report.h>
%include "Hermite.h"
%include "NonbondedMultiGridSystemForce.h"
//%include "NonbondedPMEwaldSystemForce.h"

%extend ProtoMol::NonbondedMultiGridSystemForce {
   Force* makeNew(Real s, unsigned int gridsize, int levels, unsigned int order, unsigned int ratio) {
      std::vector<Value> v;
      v.push_back(Value(s));
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(levels));
      v.push_back(Value(order));
      v.push_back(Value(ratio));
      return self->make(v);
   }
   Force* makeNew(Real s, Real h, Real o, int levels, unsigned int order, unsigned int ratio) {
      std::vector<Value> v;
      v.push_back(Value(s));
      v.push_back(Value(Vector3D(h,h,h)));
      v.push_back(Value(Vector3D(o,o,o)));
      v.push_back(Value(levels));
      v.push_back(Value(order));
      v.push_back(Value(ratio));
      return self->make(v);
   }
}




%template(MG_P_TTT) ProtoMol::NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true,true>;
%template(MG_P_FTT) ProtoMol::NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,true>;
%template(MG_P_FTF) ProtoMol::NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,false>;
%template(MG_P_FFT) ProtoMol::NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,false,true>;

%template(MG_V_TTT) ProtoMol::NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true,true>;
%template(MG_V_FTT) ProtoMol::NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,true>;
%template(MG_V_FTF) ProtoMol::NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,false>;
%template(MG_V_FFT) ProtoMol::NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,false,true>;




