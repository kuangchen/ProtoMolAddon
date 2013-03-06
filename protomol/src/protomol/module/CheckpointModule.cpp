#include "CheckpointModule.h"

#include "IOModule.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/io/CheckpointConfigReader.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/factory/OutputFactory.h>

using namespace std;
using namespace ProtoMol;

void CheckpointModule::init(ProtoMolApp *app) {
  OutputFactory &f = app->outputFactory;
  f.registerExemplar(new OutputCheckpoint());
}

void CheckpointModule::configure(ProtoMolApp *app) {
  Configuration &config = app->config;

  if (!(enabled = config.valid("Checkpoint"))) return;

  config["CheckpointPosBase"] = WithoutExt(config[InputPositions::keyword]);
  
  if (config.valid(InputVelocities::keyword))
    config["CheckpointVelBase"] = WithoutExt(config[InputVelocities::keyword]);
  else config["CheckpointVelBase"] = config["CheckpointPosBase"];
}

void CheckpointModule::read(ProtoMolApp *app) {
  if (!enabled) return;

  Configuration &config = app->config;
  
  CheckpointConfigReader confReader;
  if (confReader.open(config["Checkpoint"], ios::in))
    confReader.readBase(config, Random::Instance());
}

void CheckpointModule::postBuild(ProtoMolApp *app) {
  if (!enabled) return;

  Configuration &config = app->config;
  
  // Load integrator data
  CheckpointConfigReader confReader;
  if (confReader.open(config["Checkpoint"], ios::in))
    confReader.readIntegrator(app->integrator);
}

string CheckpointModule::WithoutExt(const string &path){
  return path.substr(0, path.rfind(".") + 1);
}
