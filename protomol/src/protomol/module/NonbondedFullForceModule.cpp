#include <protomol/module/NonbondedFullForceModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/force/OneAtomPairFull.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/LennardJonesForce.h>

#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>

#include <protomol/force/nonbonded/NonbondedFullSystemForce.h>

#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/topology/CellListEnumeratorStandard.h>

using namespace std;
using namespace ProtoMol;

void NonbondedFullForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef C1SwitchingFunction C1;
  typedef C2SwitchingFunction C2;
  typedef CnSwitchingFunction Cn;
  typedef CmpCnCnSwitchingFunction CmpCnCn;
#define FullSystem NonbondedFullSystemForce

  if (equalNocase(boundConds, PBC::keyword)) {
    f.reg(new FullSystem<OneAtomPairFull<PBC, C1, CoulombForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, C1, LennardJonesForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, C2, CoulombForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, C2, LennardJonesForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, Cn, CoulombForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, Cn, LennardJonesForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, CmpCnCn, CoulombForce> >());
    f.reg(new FullSystem<OneAtomPairFull<PBC, CmpCnCn, LennardJonesForce> >());

  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
  }
}
