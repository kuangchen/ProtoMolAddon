%module OpenMMIntegrator
%{
#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/integrator/StandardIntegrator.h>
#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/MTSIntegrator.h>
#include "ndarrayobject.h"
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%feature ("dynamic_cast");
%include "std_string.i"
%include <protomol/integrator/Integrator.h>
%include <protomol/integrator/StandardIntegrator.h>
%include <protomol/integrator/STSIntegrator.h>
%include <protomol/integrator/MTSIntegrator.h>
%include <protomol/integrator/openMM/OpenMMIntegrator.h>

%extend ProtoMol::OpenMMIntegrator {
ProtoMolApp* appInit(GenericTopology* topo,
                     Vector3DBlock& positions,
                      Vector3DBlock& velocities,
                      ScalarStructure energies) {
   import_array1(NULL);
   ProtoMolApp* app = new ProtoMolApp();
   app->topology = topo;
   app->positions.vec = positions.vec;
   app->positions.c = positions.c;
   app->velocities.vec = velocities.vec;
   app->velocities.c = velocities.c;
   app->energies = energies;
   self->initialize(app);
   app->integrator = self;
   app->outputCache.initialize(app);
   return app;
}

OpenMMIntegrator* setParams(float ts, float temperature, float gamma, int seed, 
                            bool hbf, bool haf, bool rbdf, bool pt, bool nf,
			    int it, float gbsae, float gbsaf, int com, ForceGroup* fg) {
   std::vector<Value> v;
   v.push_back(Value(ts));
   v.push_back(Value(temperature));
   v.push_back(Value(gamma));
   v.push_back(Value(seed));
   v.push_back(Value(hbf));
   v.push_back(Value(haf));
   v.push_back(Value(rbdf));
   v.push_back(Value(pt));
   v.push_back(Value(nf));
   v.push_back(Value(it));
   v.push_back(Value(gbsae));
   v.push_back(Value(gbsaf));
   v.push_back(Value(com));
   return (OpenMMIntegrator*)(self->make(v, fg)); 
}
};
