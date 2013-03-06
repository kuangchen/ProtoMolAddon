#include "VesselForcesModule.h"

#include <protomol/force/bonded/DihedralSystemForce.h>
#include "../force/VesselForce.h"
#include "../force/SlimeForce.h"
#include "../force/LJSCEExcludeForce.h"
#include "../force/PlaneForce.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

//
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/topology/CellListEnumeratorStandard.h>
#include <protomol/force/nonbonded/NonbondedCutoffSystemForce.h>
#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>

using namespace std;
using namespace ProtoMol;

void VesselForcesModule::init(ProtoMolApp *app) {}

void VesselForcesModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;

  string boundConds =
    app->config[InputBoundaryConditions::keyword];

  //Periodic BC?
  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    //Vessel wall
    f.registerExemplar(new VesselForce<PeriodicBoundaryConditions>());
    //Plane
    f.registerExemplar(new PlaneForce<PeriodicBoundaryConditions>());

    //slime
    f.registerExemplar(new SlimeForce<PeriodicBoundaryConditions>());

    //Inter-SCE LJ, C2 and Cn switching
    f.reg(new NonbondedCutoffSystemForce<CubicCellManager,
            OneAtomPair<PeriodicBoundaryConditions,
                C2SwitchingFunction, LJSCEExcludeForce> >());
    f.reg(new NonbondedCutoffSystemForce<CubicCellManager,
            OneAtomPair<PeriodicBoundaryConditions,
                CnSwitchingFunction, LJSCEExcludeForce> >());


  //Vacuum?
  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {

    //Vessel wall
    f.registerExemplar(new VesselForce<VacuumBoundaryConditions>());
    //Plane
    f.registerExemplar(new PlaneForce<VacuumBoundaryConditions>());

    //slime
    f.registerExemplar(new SlimeForce<VacuumBoundaryConditions>());

    //Inter-SCE LJ, C2 and Cn switching
    f.reg(new NonbondedCutoffSystemForce<CubicCellManager,
            OneAtomPair<VacuumBoundaryConditions,
                C2SwitchingFunction, LJSCEExcludeForce> >());
    f.reg(new NonbondedCutoffSystemForce<CubicCellManager,
            OneAtomPair<VacuumBoundaryConditions,
                CnSwitchingFunction, LJSCEExcludeForce> >());


  }
  
}
