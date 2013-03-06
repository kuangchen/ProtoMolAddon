#include <protomol/module/MainModule.h>

#include <protomol/module/IOModule.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/config/Configuration.h>
#include <protomol/base/Report.h>
#include <protomol/base/Exception.h>
#include <protomol/parallel/Parallel.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;


defineInputValue(InputSeed, "seed");
defineInputValue(InputRandomType, "randomtype");
defineInputValue(InputFirststep, "firststep");
defineInputValue(InputRealfirststep, "realfirststep");
defineInputValue(InputNumsteps, "numsteps");
defineInputValueAndText(InputDebug, "debug", "report level, suppresses all "
                        "output with higher output level");
defineInputValue(InputIntegrator, "integrator");
defineInputValue(InputReducedImage, "reducedImage");
defineInputValue(InputTemperature, "temperature");
defineInputValue(InputDoSCPISM, "doscpism");
defineInputValue(InputDoGBSAObc, "doGBSAObc");
defineInputValueAndText(InputVirialCalc, "virialCalc",
                        "Required for constant pressure simulations.");
defineInputValueAndText(InputMolVirialCalc, "molVirialCalc",
                        "Required for constant pressure simulations.");
// TODO These should be in output module
defineInputValue(InputOutputfreq,"outputfreq");
defineInputValue(InputOutput,"output");
defineInputValueAndText(InputMinimalImage, "minimalImage",
                        "global default flag whether the coordinates should be "
                        "transformed to minimal image or not");
defineInputValue(InputDebugLimit, "debugstart");


void MainModule::init(ProtoMolApp *app) {
  Configuration *config = &app->config;

  InputSeed::registerConfiguration(config, getTimerSeed());
  InputRandomType::registerConfiguration(config, 0);
  InputFirststep::registerConfiguration(config, 0);
  InputRealfirststep::registerConfiguration(config, 0);
  InputNumsteps::registerConfiguration(config);
  InputDebug::registerConfiguration(config, 1);
  InputIntegrator::registerConfiguration(config);
  InputReducedImage::registerConfiguration(config);
  InputTemperature::registerConfiguration(config);
  InputDoSCPISM::registerConfiguration(config, 0);
  InputVirialCalc::registerConfiguration(config, false);
  InputMolVirialCalc::registerConfiguration(config, false);
  InputOutput::registerConfiguration(&app->config, true);
  InputOutputfreq::registerConfiguration(&app->config, 1);
  InputMinimalImage::registerConfiguration(&app->config, false);
  InputDoGBSAObc::registerConfiguration(config, 0);
  InputDebugLimit::registerConfiguration(config, 0);
}


void MainModule::configure(ProtoMolApp *app) {
  Configuration &config = app->config;

  //  Set report level
  report << reportlevel((int)config[InputDebug::keyword],
                            (int)config[InputDebugLimit::keyword]);

  // Set random seed
  int seed = config[InputSeed::keyword];
  Parallel::bcast(seed);
  config[InputSeed::keyword] = seed;

  int randomtype;

  if (config.valid("Checkpoint")) {
    randomtype = 1;
  }else{
    randomtype = config[InputRandomType::keyword];
  }

  Parallel::bcast(randomtype);
  config[InputRandomType::keyword] = randomtype;
  randomNumber(seed, randomtype);


  // Check if configuration is complete
  if (config.hasUndefinedKeywords()) {
    report << debug(2) << "Undefined Keyword(s):" << endl
           << config.printUndefinedKeywords() << endr;
  }

  if (!config[InputFirststep::keyword].valid())
    THROW("Firststep undefined.");

  if (!config[InputNumsteps::keyword].valid())
    THROW("Numsteps undefined.");
}


void MainModule::postBuild(ProtoMolApp *app) {
  // Reduce image
  app->topology->minimalMolecularDistances =
    app->topology->checkMoleculePairDistances(app->positions);

  if ((bool)app->config[InputReducedImage::keyword] &&
      !app->topology->minimalMolecularDistances) {
    Vector3DBlock tmp(app->positions);

    app->topology->minimalImage(tmp);

    if (app->topology->checkMoleculePairDistances(tmp)) {
      app->positions = tmp;
      report << plain << "Fixed minimal molecule distances." << endr;
      app->topology->minimalMolecularDistances = true;

    } else {
      report << plain << "Could not fixed minimal molecule distances." << endr;
      app->topology->minimalMolecularDistances = false;
    }
  }

   // Fix velocities
  if (!app->config.valid(InputVelocities::keyword)) {
    randomVelocity(app->config[InputTemperature::keyword],
                   app->topology, &app->velocities,
                   app->config[InputSeed::keyword]);

    report << plain << "Random temperature : "
           << temperature(app->topology, &app->velocities) << "K" << endr;
  }
}
