#include <protomol/module/NonbondedSimpleFullForceModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomPairNoExclusion.h>
#include <protomol/force/OneAtomPairTwo.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>

#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/ComplementSwitchingFunction.h>
#include <protomol/switch/UniversalSwitchingFunction.h>

#include <protomol/force/nonbonded/NonbondedSimpleFullSystemForce.h>

#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/topology/CellListEnumeratorStandard.h>

//GB
#include <protomol/force/GB/GBBornBurialTerm.h>
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/force/GB/GBForce.h>
#include <protomol/force/GB/GBACEForce.h>


using namespace std;
using namespace ProtoMol;

void NonbondedSimpleFullForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef C1SwitchingFunction C1;
  typedef C2SwitchingFunction C2;
  typedef CnSwitchingFunction Cn;
  typedef CmpCnCnSwitchingFunction CmpCnCn;
  typedef UniversalSwitchingFunction Universal;
#define Complement ComplementSwitchingFunction 
#define SimpleFullSystem NonbondedSimpleFullSystemForce

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {

    // NonbondedSimpleFullSystemForce CoulombForce
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, C1, CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, C2, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<C1>,
          CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<C2>,
          CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<Cn>,
          CoulombForce> >());
    
    // NonbondedSimpleFullSystemForce LennardJonesForce
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Universal,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, C2, LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Cn, LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<C1>,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<C2>,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Complement<Cn>,
          LennardJonesForce> >());
    
    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Universal,
          LennardJonesForce, Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Complement<C2>,
          LennardJonesForce, Complement<C1>, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Complement<Cn>,
          LennardJonesForce, Complement<C1>, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, C2, LennardJonesForce,
          Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, C2, LennardJonesForce,
          C1, CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, C2, LennardJonesForce,
          C2, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Cn, LennardJonesForce,
          Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Cn, LennardJonesForce,
          Cn, CoulombForce> >());

    // SCPISM
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, C1, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, C2, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, Cn, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<PBC, CmpCnCn,
          CoulombSCPISMForce> >());

    // OneAtomPairTwo
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<PBC, CmpCnCn, LennardJonesForce,
          C1, CoulombSCPISMForce> >());

  } else if (equalNocase(boundConds,  VacuumBoundaryConditions::keyword)) {

    // NonbondedSimpleFullSystemForce CoulombForce
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C1, CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C2, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C1>,
          CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C2>,
          CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<Cn>,
          CoulombForce> >());

    // NonbondedSimpleFullSystemForce CoulombForceDiElec
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Universal,
          CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C1, CoulombForceDiElec> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C2, CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C1>,
          CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C2>,
          CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<Cn>,
          CoulombForceDiElec> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Universal,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C2, LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Cn, LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C1>,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<C2>,
          LennardJonesForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Complement<Cn>,
          LennardJonesForce> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Universal,
          LennardJonesForce, Universal, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2,
          LennardJonesForce, C1, CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2,
          LennardJonesForce, C2, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Complement<C2>,
          LennardJonesForce, Complement<C1>, CoulombForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Complement<Cn>,
          LennardJonesForce, Complement<C1>, CoulombForce> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForceDiElec
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Universal,
          LennardJonesForce, Universal, CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2,
          LennardJonesForce, C1, CoulombForceDiElec> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2,
          LennardJonesForce, C2, CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Complement<C2>,
          LennardJonesForce, Complement<C1>, CoulombForceDiElec> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Complement<Cn>,
          LennardJonesForce, Complement<C1>, CoulombForceDiElec> >());

    // SCPISM
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C1, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, C2, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, Cn, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPair<VBC, CmpCnCn,
          CoulombSCPISMForce> >());

    // GB
    f.reg(new SimpleFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBBornBurialTerm> >());
    f.reg(new SimpleFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBBornRadii> >());
    f.reg(new SimpleFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBACEForce> >());


    // OneAtomPairTwo
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, CmpCnCn, LennardJonesForce,
          C1, CoulombSCPISMForce> >());
    f.reg(new SimpleFullSystem<OneAtomPairTwo<VBC, Cn, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
  }
}
