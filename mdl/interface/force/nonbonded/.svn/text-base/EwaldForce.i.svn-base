%module EwaldForce
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
#include "NonbondedFullEwaldSystemForce.h"

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
%include "NonbondedFullEwaldSystemForce.h"

%extend ProtoMol::NonbondedFullEwaldSystemForce {
   Force* makeNew(Real alpha, Real accuracy) {
      std::vector<Value> v;
      v.push_back(Value(alpha));
      v.push_back(Value(accuracy));
      return self->make(v);
   }
   Force* makeNew(Real alpha, Real accuracy, Real expansionFactor) {
      std::vector<Value> v;
      v.push_back(Value(alpha));
      v.push_back(Value(accuracy));
      v.push_back(Value(expansionFactor));
      return self->make(v);
   }
}
%template(EWALD_P_TTT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,true,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_TTF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,true,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_TFF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,false,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_TFT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,false,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_FTF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,false,true,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_FTT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,false,true,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_FFT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,false,false,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_P_TTT_C1) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,true,true,ProtoMol::C1SwitchingFunction>;
%template(EWALD_P_TFT_C1) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,true,false,true,ProtoMol::C1SwitchingFunction>;
%template(EWALD_P_FTF_C1) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CubicCellManager,false,true,false,ProtoMol::C1SwitchingFunction>;



%template(EWALD_V_TTT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,true,true,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_TTF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,true,true,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_TFF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,true,false,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_TFT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,true,false,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_FTF_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,false,true,false,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_FTT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,false,true,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_FFT_CUTOFF) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,false,false,true,ProtoMol::CutoffSwitchingFunction>;
%template(EWALD_V_TTT_C1) ProtoMol::NonbondedFullEwaldSystemForce<ProtoMol::VacuumBoundaryConditions,ProtoMol::CubicCellManager,true,true,true,ProtoMol::C1SwitchingFunction>;






