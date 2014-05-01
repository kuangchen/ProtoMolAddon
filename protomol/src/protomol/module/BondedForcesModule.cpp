#include <protomol/module/BondedForcesModule.h>

#include <protomol/force/bonded/RBDihedralSystemForce.h>
#include <protomol/force/bonded/DihedralSystemForce.h>
#include <protomol/force/bonded/BondSystemForce.h>
#include <protomol/force/bonded/AngleSystemForce.h>
#include <protomol/force/bonded/ImproperSystemForce.h>
#include <protomol/force/bonded/HarmDihedralSystemForce.h>

#include <protomol/addon/ion_trap/HarmonicTrapForce.h>
#include <protomol/addon/ion_trap/RadialQuenching.h>
#include <protomol/addon/ion_trap/LQT.h>
#include <protomol/addon/template/ForceTemplate.h>
#include <protomol/addon/damping/SimpleDamping.h>
#include <protomol/addon/damping/LaserCoolingDamping.h>

#include <protomol/addon/stray_field/StrayFieldForce.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/addon/ToFForce.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::StrayField;


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
    f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::LQT>());
    f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::RadialQuenching>());
    f.registerExemplar(new IonTrap::HarmonicTrapForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new StrayFieldForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new ToFForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::SimpleDamping>());
    f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::LaserCoolingDamping>());

  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    f.registerExemplar(new RBDihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new DihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new BondSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new AngleSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new ImproperSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new HarmDihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::LQT>());
    f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::RadialQuenching>());
    f.registerExemplar(new IonTrap::HarmonicTrapForce<VacuumBoundaryConditions>());
    f.registerExemplar(new ToFForce<VacuumBoundaryConditions>());
    f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::SimpleDamping>());
    f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::LaserCoolingDamping>());
    f.registerExemplar(new StrayFieldForce<VacuumBoundaryConditions>());
  }
}
