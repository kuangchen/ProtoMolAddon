#include <protomol/module/ModifierModule.h>

#include <protomol/modifier/ModifierIncrementTimestep.h>
#include <protomol/modifier/ModifierRattle.h>
#include <protomol/modifier/ModifierShake.h>
#include <protomol/modifier/ModifierShadow.h>
#include <protomol/modifier/ModifierRemoveAngularMomentum.h>
#include <protomol/modifier/ModifierRemoveLinearMomentum.h>

#include <protomol/base/Report.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/topology/TopologyUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

defineInputValueWithAliasesAndText
(InputRemoveLinearMomentum, "removeLinearMomentum", ("comMotion"),
  "removes linear momentum, where -1 for never, 0 at initialization or at STS "
  "frequency <n>");
defineInputValueWithAliasesAndText
(InputRemoveAngularMomentum, "removeAngularMomentum", ("angularMomentum"),
  "removes angular momentum, where -1 for never, 0 at initialization or at STS "
  "frequency <n>");
defineInputValue(InputShake, "shake");
defineInputValue(InputShakeEpsilon, "shakeEpsilon");
defineInputValue(InputShakeMaxIter, "shakeMaxIter");
defineInputValue(InputShakeAll, "shakeAll");
defineInputValue(InputRattle, "rattle");
defineInputValue(InputRattleEpsilon, "rattleEpsilon");
defineInputValue(InputRattleMaxIter, "rattleMaxIter");
defineInputValue(InputRattleAll, "rattleAll");
defineInputValue(InputShadow, "shadow");
defineInputValue(InputShadowOrder, "shadoworder");
defineInputValue(InputShadowFreq, "shadowfreq");

void ModifierModule::init(ProtoMolApp *app) {
  InputRemoveLinearMomentum::registerConfiguration(&app->config, -1);
  InputRemoveAngularMomentum::registerConfiguration(&app->config, -1);
  InputShake::registerConfiguration(&app->config, false);
  InputShakeEpsilon::registerConfiguration(&app->config, 1e-5);
  InputShakeMaxIter::registerConfiguration(&app->config, 30);
  InputShakeAll::registerConfiguration(&app->config, true);
  InputRattle::registerConfiguration(&app->config, false);
  InputRattleEpsilon::registerConfiguration(&app->config, 1e-5);
  InputRattleMaxIter::registerConfiguration(&app->config, 30);
  InputRattleAll::registerConfiguration(&app->config, true);
  InputShadow::registerConfiguration(&app->config, false);
  InputShadowOrder::registerConfiguration(&app->config, 2);
  InputShadowFreq::registerConfiguration(&app->config, 1);
}

void ModifierModule::postBuild(ProtoMolApp *app) {
    // Remove Linear Momentum
  if ((int)app->config[InputRemoveLinearMomentum::keyword] >= 0)
    report
      << plain << "Removed linear momentum : "
      << removeLinearMomentum(&app->velocities, app->topology) *
      Constant::INV_TIMEFACTOR << endr;

  else
    report
      << plain << "Linear momentum : "
      << linearMomentum(&app->velocities, app->topology) *
      Constant::INV_TIMEFACTOR << endr;

  // Remove Angular Momentum
  if ((int)app->config[InputRemoveAngularMomentum::keyword] >= 0)
    report
      << plain << "Removed angular momentum : "
      << removeAngularMomentum(&app->positions, &app->velocities,
                               app->topology) *
      Constant::INV_TIMEFACTOR << endr;
  
  else
    report
      << plain << "Angular momentum : "
      << angularMomentum(&app->positions, &app->velocities, app->topology) *
      Constant::INV_TIMEFACTOR << endr;
}

void ModifierModule::addModifiers(ProtoMolApp *app) {
  Modifier *modifier;

  // Remove Linear Momentum
  int removeLinearMomentum = app->config[InputRemoveLinearMomentum::keyword];
  if (removeLinearMomentum > 0) {
    modifier = new ModifierRemoveLinearMomentum(removeLinearMomentum);
    app->integrator->bottom()->adoptPreDriftOrNextModifier(modifier);

    report << plain << "Removing linear momentum with STS frequency "
           << removeLinearMomentum << "." << endr;
  }

  // Remove Angular Momentum
  int removeAngularMomentum = app->config[InputRemoveAngularMomentum::keyword];
  if (removeAngularMomentum > 0) {
    modifier = new ModifierRemoveAngularMomentum(removeAngularMomentum);
    app->integrator->bottom()->adoptPreDriftOrNextModifier(modifier);

    report << plain << "Removing angular momentum with STS frequency "
           << removeAngularMomentum << "." << endr;
  }

  // Shake
  bool shake = app->config[InputShake::keyword];
  Real shakeEpsilon = app->config[InputShakeEpsilon::keyword];
  int shakeMaxIter = app->config[InputShakeMaxIter::keyword];
  bool shakeAll = app->config[InputShakeAll::keyword];
  if (shake && shakeEpsilon > 0.0 && shakeMaxIter > 0) {
    modifier = new ModifierShake(shakeEpsilon, shakeMaxIter, shakeAll);
    app->integrator->bottom()->adoptPostDriftOrNextModifier(modifier);

    report << plain << "Shake with epsilon " << shakeEpsilon << ", max "
           << shakeMaxIter << " iteration(s)." << endr;    
  }   

  // Rattle
  bool rattle = app->config[InputRattle::keyword];
  Real rattleEpsilon = app->config[InputRattleEpsilon::keyword];
  int rattleMaxIter = app->config[InputRattleMaxIter::keyword];
  bool rattleAll = app->config[InputRattleAll::keyword];
  if (rattle && rattleEpsilon > 0.0 && rattleMaxIter > 0) {
    modifier = new ModifierRattle(rattleEpsilon, rattleMaxIter, rattleAll);
    app->integrator->bottom()->adoptPostStepModifier(modifier);

    report << plain << "Rattle with epsilon " << rattleEpsilon <<", max "
           << rattleMaxIter << " iteration(s)." << endr;
  }

  // Shadow
  bool shadow = app->config[InputShadow::keyword];
  int shadowOrder = app->config[InputShadowOrder::keyword];
  int shadowFreq = app->config[InputShadowFreq::keyword];
  
  if (shadow) {
    modifier = new ModifierShadow(shadowOrder, shadowFreq, app->integrator);
    app->integrator->bottom()->adoptPostStepModifier(modifier);
    report << plain << "Calculating shadow energy with order 2k = "
             << shadowOrder << ", and frequency = "
             << shadowFreq << "."
             << endr;
  }
}
