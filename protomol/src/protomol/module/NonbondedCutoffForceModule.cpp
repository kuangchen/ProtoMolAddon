#include <protomol/module/NonbondedCutoffForceModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomPairTwo.h>
#include <protomol/force/OneAtomPairThree.h>
#include <protomol/force/OneAtomPairNoExclusion.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/born/BornRadii.h>
#include <protomol/force/born/BornSelfForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffSystemForce.h>
#include <protomol/force/table/LennardJonesTableForce.h>

#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/UniversalSwitchingFunction.h>
#include <protomol/switch/CutoffSwitchingFunction.h>

#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/topology/CellListEnumeratorStandard.h>

//GB
#include <protomol/force/GB/GBForce.h>
#include <protomol/force/GB/GBACEForce.h>

using namespace std;
using namespace ProtoMol;

void NonbondedCutoffForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  // To make this a bit more readable
  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef CubicCellManager CCM;
  typedef C1SwitchingFunction C1;
  typedef C2SwitchingFunction C2;
  typedef CnSwitchingFunction Cn;
  typedef CmpCnCnSwitchingFunction CmpCnCn;
  typedef UniversalSwitchingFunction Universal;
  typedef CutoffSwitchingFunction Cutoff;
#define CutoffSystem NonbondedCutoffSystemForce

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    // NonbondedCutoffSystemForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C1, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C2, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, CmpCnCn, CoulombForce> >());
    
    // NonbondedCutoffSystemForce LennardJonesForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C1, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C2, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Cn, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, CmpCnCn,
          LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<C2, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<Cn, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<CmpCnCn, 7 ,Real> > >());
    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce,
          C1, CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce,
          C2, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce,
          Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, Cn, LennardJonesForce,
          Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, CmpCnCn, LennardJonesForce,
          CmpCnCn, CoulombForce> >());


    // SCPISM stuff
    // NonbondedCutoffSystemForce CoulombSCPISMForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C1, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C2, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Cn, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, CmpCnCn,
          CoulombSCPISMForce> >());

    // OneAtomPairTwo
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, CmpCnCn, LennardJonesForce,
          C1, CoulombSCPISMForce> >());


  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    // NonbondedCutoffSystemForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn, CoulombForce> >());

    // NonbondedCutoffSystemForce LennardJonesForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn,
          LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<C2, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<Cn, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<CmpCnCn, 7 ,Real> > >());
    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C1,
          CoulombForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:

    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C2,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, Cn,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, Cn, LennardJonesForce, Cn,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, CmpCnCn, LennardJonesForce,
          C1, CoulombForce> >());

    // SCPISM stuff
    // Born radius
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cutoff, BornRadii> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cutoff, BornSelfForce> >());

    // CoulombSCPISMForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn,
          CoulombSCPISMForce> >());
    
    // GB
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, Universal, GBForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, Universal, GBACEForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, C2, GBForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, C2, GBACEForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, Cn, GBForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairNoExclusion<VBC, Cn, GBACEForce> >());
    
    // OneAtomPairThree (forces): LennardJonesForce CoulombSCPISMForce BornRadii

    f.reg(new CutoffSystem<CCM, OneAtomPairThree<VBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce, Cutoff, BornRadii > >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new CutoffSystem<CCM, OneAtomPairThree<VBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce, Cutoff, BornRadii > >());
    f.reg(new CutoffSystem<CCM, OneAtomPairThree<VBC, C2, LennardJonesForce, Cn,
          CoulombSCPISMForce, Cutoff, BornRadii > >());
    f.reg(new CutoffSystem<CCM, OneAtomPairThree<VBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce, Cutoff, BornRadii > >());


    // OneAtomPairTwo
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce> >());
    // this option should be used with switchon = 0 for Coulomb
    // C2 continuity needed when computing Hessians:
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, CmpCnCn, LennardJonesForce,
          C1, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, Cn, LennardJonesForce, C2,
          CoulombSCPISMForce> >());


    // End of SCPISM



  }
}
