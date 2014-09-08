#include <protomol/addon/AddonModule.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/addon/LeapfrogBufferGasIntegrator2.h>
#include <protomol/addon/reaction/ReactionIntegrator.h>
#include <protomol/addon/reaction/Replacement.h>
#include <protomol/addon/reaction/Ionization.h>
#include <protomol/addon/reaction/EscapeDetection.h>

#include <protomol/addon/buffergas/LeapfrogBufferGasIntegrator.h>
#include <protomol/addon/buffergas/IsotropicCollision.h>
#include <protomol/addon/damping/SimpleDamping.h>

#include <protomol/addon/ion_trap/HarmonicTrapForce.h>
#include <protomol/addon/ion_trap/RadialQuenching.h>
#include <protomol/addon/ion_trap/LQT.h>
#include <protomol/addon/template/ForceTemplate.h>
#include <protomol/addon/damping/SimpleDamping.h>
#include <protomol/addon/damping/LaserCoolingDamping.h>
#include <protomol/addon/stray_field/StrayField.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/addon/template/ForceTemplate.h>

#include <protomol/addon/tof/OutputCEMRecorder.h>

#include <protomol/addon/snapshot/OutputSnapshot.h>
#include <protomol/addon/snapshot/HDF5Storage.h>


using namespace ProtoMol;

namespace ProtoMolAddon {

  void AddonModule::init(ProtoMol::ProtoMolApp *app) {
    app->integratorFactory.registerExemplar(new LeapfrogBufferGasIntegrator2());
    app->integratorFactory.registerExemplar(new Reaction::ReactionIntegrator<Reaction::Replacement>);
    app->integratorFactory.registerExemplar(new Reaction::ReactionIntegrator<Reaction::Ionization>);
    app->integratorFactory.registerExemplar(new Reaction::ReactionIntegrator<Reaction::EscapeDetection>);

    app->integratorFactory.registerExemplar(new BufferGas::LeapfrogBufferGasIntegrator<BufferGas::IsotropicCollision>);

    OutputFactory &f = app->outputFactory;  
    f.registerExemplar(new ToF::OutputCEMRecorder()); 
    f.registerExemplar(new Snapshot::OutputSnapshot<Snapshot::HDF5Storage>()); 
  }




  void AddonModule::registerForces(ProtoMol::ProtoMolApp *app) {
    ForceFactory &f = app->forceFactory;

    string boundConds =
      app->config[InputBoundaryConditions::keyword];

    if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::SimpleDamping>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, Damping::LaserCoolingDamping>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::LQT>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, IonTrap::RadialQuenching>());
      f.registerExemplar(new Template::GenericForce<PeriodicBoundaryConditions, StrayField::StrayField>());
      f.registerExemplar(new IonTrap::HarmonicTrapForce<PeriodicBoundaryConditions>());


    } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::SimpleDamping>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, Damping::LaserCoolingDamping>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::LQT>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, IonTrap::RadialQuenching>());
      f.registerExemplar(new Template::GenericForce<VacuumBoundaryConditions, StrayField::StrayField>());
      f.registerExemplar(new IonTrap::HarmonicTrapForce<VacuumBoundaryConditions>());
    }
  }

}


