#include <protomol/module/BondedForcesModule.h>

#include <protomol/force/bonded/RBDihedralSystemForce.h>
#include <protomol/force/bonded/DihedralSystemForce.h>
#include <protomol/force/bonded/BondSystemForce.h>
#include <protomol/force/bonded/AngleSystemForce.h>
#include <protomol/force/bonded/ImproperSystemForce.h>
#include <protomol/force/bonded/HarmDihedralSystemForce.h>

#include <protomol/addon/ion_trap/HarmonicTrapForce.h>
#include <protomol/addon/LQTForce.h>
//#include <protomol/addon/LaserCoolingForce.h>
#include <protomol/addon/DampingForce.h>
#include <protomol/addon/StrayFieldForce.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/addon/ToFForce.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::IonTrap;

void BondedForcesModule::init(ProtoMolApp *app) {}

void BondedForcesModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;

  string boundConds =
    app->config[InputBoundaryConditions::keyword];

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    f.registerExemplar(new RBDihedralSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new DihedralSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new BondSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new AngleSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new ImproperSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new HarmDihedralSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new LQTForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new HarmonicTrapForce<PeriodicBoundaryConditions>());
//f.registerExemplar(new LaserCoolingForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new DampingForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new StrayFieldForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new ToFForce<PeriodicBoundaryConditions>());

  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    f.registerExemplar(new RBDihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new DihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new BondSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new AngleSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new ImproperSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new HarmDihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new LQTForce<VacuumBoundaryConditions>());
    f.registerExemplar(new HarmonicTrapForce<VacuumBoundaryConditions>());
    f.registerExemplar(new DampingForce<VacuumBoundaryConditions>());
    f.registerExemplar(new ToFForce<VacuumBoundaryConditions>());
//    f.registerExemplar(new LaserCoolingForce<VacuumBoundaryConditions>());
    f.registerExemplar(new StrayFieldForce<VacuumBoundaryConditions>());
  }
}
