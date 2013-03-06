%module PMEForce
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
#include "BSpline.h"
#include "Hermite.h"
#include "NonbondedPMEwaldSystemForce.h"

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
%include "BSpline.h"
%include "Hermite.h"
%include "NonbondedPMEwaldSystemForce.h"

%extend ProtoMol::NonbondedPMEwaldSystemForce {
   Force* makeNew(unsigned int gridsize, unsigned int order, Real rc, Real alpha, Real accuracy) {
      std::vector<Value> v;
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(rc));
      v.push_back(Value(order));
      v.push_back(Value(accuracy));
      v.push_back(Value(alpha));
      return self->make(v);
   }
   Force* makeNew(unsigned int gridsize, unsigned int order, Real rc, Real alpha, Real accuracy, Real expansionFactor) {
      std::vector<Value> v;
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(gridsize));
      v.push_back(Value(rc));
      v.push_back(Value(order));
      v.push_back(Value(accuracy));
      v.push_back(Value(alpha));
      v.push_back(Value(expansionFactor));
      return self->make(v);
   }
}



%template(PME_P_TTT_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>;
%template(PME_P_TFF_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>;
%template(PME_P_TFT_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>;
%template(PME_P_FTF_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>;
%template(PME_P_FTT_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>;
%template(PME_P_FFT_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>;
%template(PME_P_TTF_B) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>;
%template(PME_P_TTT_B_C1) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>;

%template(PME_P_TTT_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>;
%template(PME_P_TFF_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>;
%template(PME_P_TFT_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>;
%template(PME_P_FTF_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>;
%template(PME_P_FTT_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>;
%template(PME_P_FFT_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>;
%template(PME_P_TTF_H) ProtoMol::NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>;


%template(PME_V_TTT_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>;
%template(PME_V_TFF_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>;
%template(PME_V_TFT_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>;
%template(PME_V_FTF_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>;
%template(PME_V_FTT_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>;
%template(PME_V_FFT_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>;
%template(PME_V_TTF_B) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>;
%template(PME_V_TTT_B_C1) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>;

%template(PME_V_TTT_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>;
%template(PME_V_TFF_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>;
%template(PME_V_TFT_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>;
%template(PME_V_FTF_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>;
%template(PME_V_FTT_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>;
%template(PME_V_FFT_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>;
%template(PME_V_TTF_H) ProtoMol::NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>;



