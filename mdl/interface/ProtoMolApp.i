%module ProtoMolApp
%{
#include <protomol/ProtoMolApp.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>
#include <protomol/base/Report.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <interface/type/Vector3DBlock.i>
%include <interface/type/ScalarStructure.i>
%include <protomol/ProtoMolApp.h>

%extend ProtoMol::ProtoMolApp {
static void turnOffHints() {
   Report::report << Report::donthint;
}
void uncache() {
   self->outputCache.uncache();
   self->config.registerKeyword("firststep", Value(0));
}
void makeApp(GenericTopology* topo,
             ProtoMol::Vector3DBlock& positions,
             ProtoMol::Vector3DBlock& velocities,
             ScalarStructure energies,
             Real timestep) {
   self->config.registerKeyword("firststep", Value(0));
   self->topology = topo;
   self->positions.vec = positions.vec;
   self->positions.c = positions.c;
   self->velocities.vec = velocities.vec;
   self->velocities.c = velocities.c;
   self->energies = energies;
   self->integrator = new LeapfrogIntegrator(timestep, NULL);
   self->outputCache.initialize(self);
}
};
