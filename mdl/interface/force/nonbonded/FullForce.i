%module NonbondedFullSystemForce
%{
#include "CoulombForce.h"
#include "LennardJonesForce.h"
#include "../topology/C1SwitchingFunction.h"
#include "../topology/C2SwitchingFunction.h"
#include "../topology/CutoffSwitchingFunction.h"
#include "../topology/ShiftSwitchingFunction.h"
#include "../topology/CnSwitchingFunction.h"
#include "../topology/CmpCnCnSwitchingFunction.h"
#include "../topology/PeriodicBoundaryConditions.h"
#include "OneAtomPairFull.h"
#include "OneAtomPairTwoFull.h"
#include "NonbondedFullSystemForce.h"

void setSwitchonFCP(ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn)
{  
   ljsf->myOneAtomPair.switchingFunction.mySwitchon = sn;
   return;
}
void setSwitchonFLP(ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >* ljsf, float sn)
{  
   ljsf->myOneAtomPair.switchingFunction.mySwitchon = sn;
   return;
}


%}

%include "std_vector.i"
%template() std::vector<ProtoMol::Vector3D >;
%include "../base/Real.h"
%include "../base/Proxy.h"
%include "../base/Vector3D.h"
%include "../base/AbstractEnumType.h"
%include "../base/simpleTypes.h"
%include "../base/ConstraintValueType.h"
%include "../base/ValueType.h"
%include "../base/Value.h"
%include "../base/Parameter.h"

%include "../base/Vector3DBlock.h"

%include "../base/MakeableDefinition.h"
%include "../base/Makeable.h"





%include "Force.h"
%include "SystemForce.h"
%include "NonbondedFullSystemForce.h"
 %template (NFSF_P_C1_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
 %template (NFSF_P_C1_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> >;
 %template (NFSF_P_C2_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >;
 %template (NFSF_P_C2_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >;
 %template (NFSF_P_CO_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;
 %template (NFSF_P_CO_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce> >;
 %template (NFSF_P_S_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ShiftSwitchingFunction,ProtoMol::CoulombForce> >;
 %template (NFSF_P_S_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ShiftSwitchingFunction,ProtoMol::LennardJonesForce> >;


 %template (NFSF_P_C1_L_S_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairTwoFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::ShiftSwitchingFunction,ProtoMol::CoulombForce> >;
 %template (NFSF_P_CO_L_C_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairTwoFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;

%template(NFSF_P_CN_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NFSF_P_CCNCN_C) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NFSF_P_CN_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NFSF_P_CCNCN_L) ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce> >;


void setSwitchonFCP(ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn);
void setSwitchonFLP(ProtoMol::NonbondedFullSystemForce<ProtoMol::OneAtomPairFull<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >* ljsf, float sn);




