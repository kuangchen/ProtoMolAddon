%module CutoffForce
%{
#include <protomol/base/Report.h>
#include <protomol/config/Value.h>
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
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/GB/GBForce.h>
#include <protomol/force/GB/GBACEForce.h>
//#include <protomol/force/coulomb/CoulombBornRadiiForce.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/table/LennardJonesTableForce.h>
#include <protomol/force/table/CoulombTableForce.h>
#include <protomol/topology/CellListEnumeratorStandard.h>
#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/force/nonbonded/NonbondedCutoffForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffSystemForce.h>
//#include <protomol/force/nonbonded/NonbondedCutoffBornForce.h>
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
%include <protomol/force/nonbonded/NonbondedCutoffForce.h>
%include <protomol/force/nonbonded/NonbondedCutoffSystemForce.h>
//%include <protomol/force/nonbonded/NonbondedCutoffBornForce.h>

%extend ProtoMol::NonbondedCutoffForce {
   Force* makeNewDiElec(Real cutoff,
			Real eps, Real d, Real s,
                        Real switchon=-1,
                        Real switchoff=-1,
                        Real order=-1) {
      std::vector<Value> v;
      if (switchon == -1) { // C1
	 v.push_back(Value(eps));
	 v.push_back(Value(d));
	 v.push_back(Value(s));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      else if (order == -1) { // C2
	 v.push_back(Value(eps));
	 v.push_back(Value(d));
	 v.push_back(Value(s));
         v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      else {
	 v.push_back(Value(eps));
	 v.push_back(Value(d));
	 v.push_back(Value(s));
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
         v.push_back(Value(cutoff));
      }
      return self->make(v);
   }


   Force* makeNewGB(Real d, Real s,
		    Real cutoff=-1,
                    Real switchon=-1,
                    Real switchoff=-1,
                    Real order=-1) {
      std::vector<Value> v;
      if (switchon == -1) { // Universal
	 v.push_back(Value(d));
	 v.push_back(Value(s));
	 v.push_back(Value(cutoff));
      }
      else if (order == -1) { // C2
	 v.push_back(Value(d));
	 v.push_back(Value(s));
         v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      else { // Cn
	 v.push_back(Value(d));
	 v.push_back(Value(s));
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
         v.push_back(Value(cutoff));
      }
      return self->make(v);
   }


   Force* makeNew(Real cutoff,
                  Real switchon=-1,
                  Real switchoff=-1,
                  Real order=-1) {
      std::vector<Value> v;
      if (switchon == -1) { // C1
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      else if (order == -1) { // C2
         v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      else { // Cn
         v.push_back(Value(order));
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));
      }
      return self->make(v);
   }

   Force* makeNewPair(std::string switch1,
		      std::string switch2,
		      Real cutoff,
                      Real switchon=-1,
		      Real switchoff=-1,
                      Real order=-1) {
      std::vector<Value> v;
      if (switch1 == "C1" && switch2 == "C1") {
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));	
      }
      else if (switch1 == "C1" && switch2 == "C2") {
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));         
      }
      else if (switch1 == "C2" && switch2 == "C1") {
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));           
      }
      else if (switch1 == "C2" && switch2 == "C2") {
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));  
      }
      else if (switch1 == "C1" && switch2 == "Cn") {
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));  
      }
      else if (switch1 == "C2" && switch2 == "Cn") {
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));	
         v.push_back(Value(cutoff));  
      }
      else if (switch1 == "Cn" && switch2 == "C1") {
         v.push_back(Value(cutoff));  
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
         v.push_back(Value(cutoff));  
      }
      else if (switch1 == "Cn" && switch2 == "C2") {
	 v.push_back(Value(switchon));
         v.push_back(Value(cutoff));
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));	
         v.push_back(Value(cutoff));  
      }
      else if (switch1 == "Cn" && switch2 == "Cn") {
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
	 v.push_back(Value(switchon));
         v.push_back(Value(switchoff));
         v.push_back(Value(cutoff));
         v.push_back(Value(order));
         v.push_back(Value(cutoff));  
      }
      return self->make(v);
}
};


//******************************
// TEMPLATES
//******************************

%template(NCF_CCM_OAPPBC_C1SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPPBC_C2SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPPBC_CSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPPBC_C1SF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NSF_CCM_OAPPBC_C2SF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPPBC_CSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPTPBC_C2SF_LJF_C1SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTPBC_C2SF_LJF_CSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;



%template(NCF_CCM_OAPVBC_C1SF_CSCPF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombSCPISMForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombSCPISMForce> > >;
%template(NCF_CCM_OAPVBC_C2SF_CSCPF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombSCPISMForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombSCPISMForce> > >;
%template(NCF_CCM_OAPVBC_CNSF_CSCPF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombSCPISMForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombSCPISMForce> > >;
%template(NCF_CCM_OAPVBC_CCNSF_CSCPF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombSCPISMForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombSCPISMForce> > >;
%template(NCF_CCM_OAPVBC_CSF_CSCPF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombSCPISMForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombSCPISMForce> > >;



%template(NCF_CCM_OAPVBC_U_GB) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBForce> > >;
%template(NCF_CCM_OAPVBC_U_GBACE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBACEForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBACEForce> > >;
%template(NCF_CCM_OAPVBC_C2_GB) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBForce> > >;
%template(NCF_CCM_OAPVBC_C2_GBACE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBACEForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBACEForce> > >;
%template(NCF_CCM_OAPVBC_CN_GB) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBForce> > >;
%template(NCF_CCM_OAPVBC_CN_GBACE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBACEForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBACEForce> > >;







//%template(NCF_CCM_OAPVBC_C1SF_CBF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombBornRadiiForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombBornRadiiForce> > >;
//%template(NCF_CCM_OAPVBC_C2SF_CBF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombBornRadiiForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombBornRadiiForce> > >;
//%template(NCF_CCM_OAPVBC_CNSF_CBF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombBornRadiiForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombBornRadiiForce> > >;
//%template(NCF_CCM_OAPVBC_CCNSF_CBF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombBornRadiiForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombBornRadiiForce> > >;
//%template(NCF_CCM_OAPVBC_CSF_CBF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombBornRadiiForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombBornRadiiForce> > >;



%template(NCF_CCM_OAPVBC_C1SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPVBC_C2SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPVBC_CSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPVBC_C1SF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NSF_CCM_OAPVBC_C2SF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPVBC_CSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPTVBC_C2SF_LJF_C1SF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTVBC_C2SF_LJF_CSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;







//******************************
// SYSTEM FORCES
//******************************


%template(NCSF_CCM_OAPPBC_C1SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPPBC_C2SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPPBC_CSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPPBC_C1SF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPPBC_C2SF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPPBC_CSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPTPBC_C2SF_LJF_C1SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTPBC_C2SF_LJF_CSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_C1SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_C2SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_CSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_C1SF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPVBC_C2SF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPVBC_CSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPTVBC_C2SF_LJF_C1SF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTVBC_C2SF_LJF_CSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;


// DIELEC
%template(NCF_CCM_OAPPBC_CNSF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec> > >;
%template(NCF_CCM_OAPVBC_C1SF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForceDiElec> > >;
%template(NCF_CCM_OAPVBC_C2SF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForceDiElec> > >;
%template(NCF_CCM_OAPVBC_CNSF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec> > >;
%template(NCF_CCM_OAPVBC_CMPCNNSF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForceDiElec> > >;
%template(NCF_CCM_OAPVBC_CSF_CFDE) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForceDiElec>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForceDiElec> > >;

// SCPISM
%template(NCSF_CCM_OAPVBC_C1SF_CSCPF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombSCPISMForce> >;
%template(NCSF_CCM_OAPVBC_C2SF_CSCPF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombSCPISMForce> >;
%template(NCSF_CCM_OAPVBC_CNSF_CSCPF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombSCPISMForce> >;
%template(NCSF_CCM_OAPVBC_CCNSF_CSCPF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombSCPISMForce> >;
%template(NCSF_CCM_OAPVBC_CSF_CSCPF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombSCPISMForce> >;


// Born
//%template(NCBF_CCM_OAPVBC_C1SF_CBF) ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombBornRadiiForce> >;
//%template(NCBF_CCM_OAPVBC_C2SF_CBF) ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombBornRadiiForce> >;
//%template(NCBF_CCM_OAPVBC_CNSF_CBF) ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombBornRadiiForce> >;
//%template(NCBF_CCM_OAPVBC_CCNSF_CBF) ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombBornRadiiForce> >;
//%template(NCBF_CCM_OAPVBC_CSF_CBF) ProtoMol::NonbondedCutoffBornForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombBornRadiiForce> >;




// ************************** CN *******************************
// SYSTEM FORCES
%template(NCF_CCM_OAPPBC_CNSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPPBC_CCNCNSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPPBC_CNSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPPBC_CCNCNSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce> > >;

%template(NCF_CCM_OAPPBC_U_LJTF_CN) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> >, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> > > >;
%template(NCF_CCM_OAPPBC_U_LJTF_CCNCN) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> >, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> > > >;

%template(NCF_CCM_OAPTPBC_CN_LJF_CN_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTPBC_CN_LJF_C_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;


%template(NCF_CCM_OAPTPBC_CCNCN_LJF_CCNCN_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTPBC_CCNCN_LJF_C_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;


%template(NCF_CCM_OAPVBC_CNSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPVBC_CCNCNSF_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPVBC_CNSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> > >;
%template(NCF_CCM_OAPVBC_CCNCNSF_LJF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce> > >;

%template(NCF_CCM_OAPVBC_U_LJTF_CN) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> >, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> > > >;
%template(NCF_CCM_OAPVBC_U_LJTF_CCNCN) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> >, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> > > >;

%template(NCF_CCM_OAPTVBC_CN_LJF_CN_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTVBC_CN_LJF_C_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;


%template(NCF_CCM_OAPTVBC_CCNCN_LJF_CCNCN_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> > >;
%template(NCF_CCM_OAPTVBC_CCNCN_LJF_C_CF) ProtoMol::NonbondedCutoffForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce>, ProtoMol::SystemForce, ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> > >;


// SYSTEM FORCES
%template(NCSF_CCM_OAPPBC_CNSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPPBC_CCNCNSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPPBC_CNSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPPBC_CCNCNSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce> >;

%template(NCSF_CCM_OAPPBC_U_LJTF_CN) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> > >;
%template(NCSF_CCM_OAPPBC_U_LJTF_CCNCN) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> > >;

%template(NCSF_CCM_OAPTPBC_CN_LJF_CN_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTPBC_CN_LJF_C_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;


%template(NCSF_CCM_OAPTPBC_CCNCN_LJF_CCNCN_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTPBC_CCNCN_LJF_C_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;



%template(NCSF_CCM_OAPVBC_CNSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_CCNCNSF_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPVBC_CNSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce> >;
%template(NCSF_CCM_OAPVBC_CCNCNSF_LJF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce> >;

%template(NCSF_CCM_OAPVBC_U_LJTF_CN) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CnSwitchingFunction,7,ProtoMol::Real> > >;
%template(NCSF_CCM_OAPVBC_U_LJTF_CCNCN) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::UniversalSwitchingFunction,ProtoMol::LennardJonesTableForce<ProtoMol::CmpCnCnSwitchingFunction,7,ProtoMol::Real> > >;

%template(NCSF_CCM_OAPTVBC_CN_LJF_CN_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTVBC_CN_LJF_C_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;


%template(NCSF_CCM_OAPTVBC_CCNCN_LJF_CCNCN_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForce> >;
%template(NCSF_CCM_OAPTVBC_CCNCN_LJF_C_CF) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >;





//******************************
// DECLARATIONS
//******************************


/*void setSwitchonCEP(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForceDiElec> >* ljsf, float sn);
void setSwitchonCCP(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn);
void setSwitchonCLP(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >* ljsf, float sn);*/


//********************************************
// VACUUM
//********************************************

/*void setSwitchonCEV(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForceDiElec> >* ljsf, float sn);
void setSwitchonCCV(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn);
void setSwitchonCLV(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce> >* ljsf, float sn);
void setSwitchonCLPPC2C1(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn1);
void setSwitchonCLPPC2Cutoff(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::PeriodicBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn1);
void setSwitchonCLVPC2C1(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn1);
void setSwitchonCLVPC2Cutoff(ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPairTwo<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::LennardJonesForce,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForce> >* ljsf, float sn1);*/




%template(NCSF_CCM_OAPPBC_CNSF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::PeriodicBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NCSF_CCM_OAPVBC_C1SF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C1SwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NCSF_CCM_OAPVBC_C2SF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::C2SwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NCSF_CCM_OAPVBC_CNSF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CnSwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NCSF_CCM_OAPVBC_CMPCNNSF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CmpCnCnSwitchingFunction,ProtoMol::CoulombForceDiElec> >;
%template(NCSF_CCM_OAPVBC_CSF_CFDE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager,ProtoMol::OneAtomPair<ProtoMol::VacuumBoundaryConditions,ProtoMol::CutoffSwitchingFunction,ProtoMol::CoulombForceDiElec> >;


# GBSA
%template(NCSF_CCM_OAPVBC_U_GB) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBForce> >;
%template(NCSF_CCM_OAPVBC_U_GBACE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::UniversalSwitchingFunction, ProtoMol::GBACEForce> >;
%template(NCSF_CCM_OAPVBC_C2_GB) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBForce> >;
%template(NCSF_CCM_OAPVBC_C2_GBACE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::C2SwitchingFunction, ProtoMol::GBACEForce> >;
%template(NCSF_CCM_OAPVBC_CN_GB) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBForce> >;
%template(NCSF_CCM_OAPVBC_CN_GBACE) ProtoMol::NonbondedCutoffSystemForce<ProtoMol::CubicCellManager, ProtoMol::OneAtomPairNoExclusion<ProtoMol::VacuumBoundaryConditions, ProtoMol::CnSwitchingFunction, ProtoMol::GBACEForce> >;


