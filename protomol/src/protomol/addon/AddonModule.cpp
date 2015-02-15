#include <protomol/addon/AddonModule.h>
#include <protomol/ProtoMolApp.h>

//#include <protomol/addon/LeapfrogBufferGasIntegrator2.h>
#include <protomol/addon/reaction/Replacement.h>

#include <protomol/addon/buffergas/LeapfrogSimpleBufferGas.h>
#include <protomol/addon/damping/SimpleDamping.h>

#include <protomol/addon/ion_trap/PatchField.h>
#include <protomol/addon/ion_trap/LQT.h>
#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <protomol/addon/template/ForceTemplate.h>
#include <protomol/addon/template/OutputTemplate.h>

#include <protomol/addon/damping/SimpleDamping.h>
#include <protomol/addon/damping/LaserCoolingDamping.h>
#include <protomol/addon/stray_field/StrayField.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/addon/template/ForceTemplate.h>
#include <protomol/addon/template/IntegratorTemplate.h>

#include <protomol/addon/reaction/Evaporation.h>
#include <protomol/addon/reaction/Replacement.h>

#include <protomol/addon/tof/CEM.h>
#include <protomol/addon/snapshot/OutputSnapshot.h>
//#include <protomol/addon/snapshot/HDF5Storage.h>

namespace ProtoMolAddon {
  using namespace ProtoMol;
  
  void AddonModule::init(ProtoMolApp *app) {
    IntegratorFactory &i = app->integratorFactory;
    
    i.registerExemplar(new Template::GenericIntegrator<Reaction::Evaporation>);
    i.registerExemplar(new Template::GenericIntegrator<Reaction::Replacement>);
    i.registerExemplar(new Template::GenericIntegrator<BufferGas::LeapfrogSimpleBufferGas>);

    OutputFactory &f = app->outputFactory;
    f.registerExemplar(new Template::GenericOutput<ToF::CEM>());
    //f.registerExemplar(new Template::GenericOutput<Snapshot::HDF5Storage>());
    
    //f.registerExemplar(new Snapshot::OutputSnapshot<Snapshot::HDF5Storage>()); 
  }

  void AddonModule::registerForces(ProtoMol::ProtoMolApp *app) {
    ForceFactory &f = app->forceFactory;

    string boundConds =
      app->config[InputBoundaryConditions::keyword];

    if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::SimpleDamping>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::LaserCoolingDamping>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::LQT>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::HarmonicTrap>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::PatchField>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, StrayField::StrayField>());

    } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::SimpleDamping>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::LaserCoolingDamping>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::LQT>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::HarmonicTrap>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::PatchField>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, StrayField::StrayField>());
    }
  }

}


