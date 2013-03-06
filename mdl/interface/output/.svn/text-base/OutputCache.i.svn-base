%module OutputCache
%{
#include <protomol/output/OutputCache.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/output/OutputCache.h>
%extend ProtoMol::OutputCache {
void appInit(GenericTopology* topo,
             Vector3DBlock positions,
             Vector3DBlock velocities,
             ScalarStructure energies) {
   ProtoMolApp* app = new ProtoMolApp();
   app->topology = topo;   app->positions = positions;   app->velocities = velocities;   app->energies = energies;   self->initialize(app);   }
};
