#include <protomol/module/NonbondedIntermittentFullForceModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/force/OneAtomPairNoExclusion.h>

#include <protomol/switch/UniversalSwitchingFunction.h>
#include <protomol/switch/CutoffSwitchingFunction.h>

#include <protomol/force/nonbonded/NonbondedIntermittentFullSystemForce.h>

#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/topology/CellListEnumeratorStandard.h>

//SCPISM
#include <protomol/force/born/BornRadii.h>
#include <protomol/force/born/BornSelfForce.h>

//GB
#include <protomol/force/GB/GBBornBurialTerm.h>
#include <protomol/force/GB/GBBornRadii.h>

using namespace std;
using namespace ProtoMol;

void NonbondedIntermittentFullForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef CutoffSwitchingFunction Cutoff;
  typedef UniversalSwitchingFunction Universal;
#define Complement ComplementSwitchingFunction 
#define IntermittentFullSystem NonbondedIntermittentFullSystemForce

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {


  } else if (equalNocase(boundConds,  VacuumBoundaryConditions::keyword)) {
    
    // SCPISM
    // Born radius
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Cutoff, BornRadii> >());
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Cutoff, BornSelfForce> >());
    
    // GB
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBBornBurialTerm> >());
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBBornRadii> >());
    
  }
}
