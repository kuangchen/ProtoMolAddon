%module SimpleFullForce
%{
#include <protomol/base/Report.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomPairTwo.h>
#include <protomol/force/OneAtomPairNoExclusion.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/switch/UniversalSwitchingFunction.h>
#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/ComplementSwitchingFunction.h>
#include <protomol/switch/CutoffSwitchingFunction.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/GB/GBForce.h>
#include <protomol/force/GB/GBACEForce.h>
#include <protomol/force/GB/GBBornBurialTerm.h>
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/force/nonbonded/NonbondedSimpleFullSystemForce.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
//%include <protomol/config/Value.h>
%include <protomol/config/Parameter.h>
%include <protomol/type/Vector3D.h>
%include <protomol/type/Vector3DBlock.h>
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
%include <protomol/force/nonbonded/NonbondedSimpleFullSystemForce.h>

%extend ProtoMol::NonbondedSimpleFullSystemForce {
   Force* makeNewDiElec(unsigned int blocksize, Real eps, Real d, Real s) {
      std::vector<Value> v;
      v.push_back(Value(eps));
      v.push_back(Value(d));
      v.push_back(Value(s));
      v.push_back(Value(blocksize));
      return self->make(v);
   }
   Force* makeNew(unsigned int blocksize=32) {
      std::vector<Value> v;
      v.push_back(Value(blocksize));
      return self->make(v);
   }
   Force* makeNewPair(unsigned int blocksize=32) {
      std::vector<Value> v;
      v.push_back(Value(blocksize));
      return self->make(v);
}
};


// NonbondedSimpleFullSystemForce CoulombForce
%template(NSFSF_P_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_P_C1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CC2_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CCO_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CutoffSwitchingFunction>,ProtoMol::CoulombForce> >;
    
    // NonbondedSimpleFullSystemForce LennardJonesForce
%template(NSFSF_P_U_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_C1_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_CC1_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_CC2_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_CCO_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CutoffSwitchingFunction>,ProtoMol::LennardJonesForce> >;
      
    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
%template(NSFSF_P_U_L_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CC2_L_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;


///VACUUM


    // NonbondedSimpleFullSystemForce CoulombForce
%template(NSFSF_V_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_V_C1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CC2_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CCO_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CutoffSwitchingFunction>,ProtoMol::CoulombForce> >;
    
    // NonbondedSimpleFullSystemForce LennardJonesForce
%template(NSFSF_V_U_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_C1_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_CC1_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_CC2_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_CCO_L) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CutoffSwitchingFunction>,ProtoMol::LennardJonesForce> >;


// DIELEC

%template(NSFSF_V_U_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NSFSF_V_C1_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NSFSF_V_CC1_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<C1SwitchingFunction>,ProtoMol::CoulombForceDiElec> >;
%template(NSFSF_V_CC2_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<C2SwitchingFunction>,ProtoMol::CoulombForceDiElec> >;
%template(NSFSF_V_CCN_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<CnSwitchingFunction>,ProtoMol::CoulombForceDiElec> >;
%template(NSFSF_V_CC_CDE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<CutoffSwitchingFunction>,ProtoMol::CoulombForceDiElec> >;

      
    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
%template(NSFSF_V_U_L_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CC2_L_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;


// Cn Switching Functions
%template(NSFSF_P_CCN_C) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CN_L) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_CCN_L) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_P_CCN_L_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_P_CN_L_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;

%template(NSFSF_V_CCN_C) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CN_L) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_CCN_L) ProtoMol::NonbondedSimpleFullSystemForce<
				     ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::LennardJonesForce> >;
%template(NSFSF_V_CCN_L_CC1_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForce> >;
%template(NSFSF_V_CN_L_U_C) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForce> >;


// DIELEC
%template(NSFSF_V_U_L_U_CDE) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::UniversalSwitchingFunction,ProtoMol::CoulombForceDiElec> >;

%template(NSFSF_V_CC2_L_CC1_CDE) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::C2SwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForceDiElec> >;


%template(NSFSF_V_CCN_L_CC1_CDE) ProtoMol::NonbondedSimpleFullSystemForce<
				   ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::ComplementSwitchingFunction<ProtoMol::CnSwitchingFunction>,ProtoMol::LennardJonesForce,ProtoMol::ComplementSwitchingFunction<ProtoMol::C1SwitchingFunction>,ProtoMol::CoulombForceDiElec> >;


%template(NSFSF_V_U_GB_GBORNBUR) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBBornBurialTerm> >;

%template(NSFSF_V_U_GB_GBORN) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBBornRadii> >;

%template(NSFSF_V_U_GB_GB) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBForce> >;

%template(NSFSF_V_U_GB_GBACE) ProtoMol::NonbondedSimpleFullSystemForce<ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBACEForce> >;
