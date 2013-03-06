#include <protomol/base/ModuleManager.h>

#include <protomol/module/MainModule.h>
#include <protomol/module/CommandLineModule.h>
#include <protomol/module/ConfigurationModule.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/module/OutputModule.h>
#include <protomol/module/BondedForcesModule.h>
#include <protomol/module/ModifierModule.h>
#include <protomol/module/IOModule.h>

#include <protomol/module/IntegratorBaseModule.h>
#include <protomol/module/NormalModeModule.h>
#include <protomol/module/LeapfrogModule.h>
#include <protomol/module/HessianIntegratorModule.h>

#include <protomol/module/NonbondedCutoffForceModule.h>
#include <protomol/module/NonbondedFullForceModule.h>
#include <protomol/module/NonbondedSimpleFullForceModule.h>

using namespace ProtoMol;

void moduleInitFunction(ModuleManager *manager) {
  manager->add(new MainModule());
  manager->add(new CommandLineModule());
  manager->add(new ConfigurationModule());
  manager->add(new TopologyModule());
  manager->add(new OutputModule());
  //manager->add(new ModifierModule());
  manager->add(new IOModule());

  // Integrators
  //manager->add(new IntegratorBaseModule());
  //manager->add(new NormalModeModule());
  manager->add(new LeapfrogModule());
  //manager->add(new HessianIntegratorModule());

  // Forces
  //manager->add(new BondedForcesModule());
  //manager->add(new NonbondedCutoffForceModule());
  //manager->add(new NonbondedFullForceModule());
  //manager->add(new NonbondedSimpleFullForceModule());
}

