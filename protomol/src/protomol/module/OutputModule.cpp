#include <protomol/module/OutputModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/config/Configuration.h>
#include <protomol/factory/OutputFactory.h>

#include <protomol/output/OutputDCDTrajectory.h>
#include <protomol/output/OutputDCDTrajectoryVel.h>
#include <protomol/output/OutputFinalPDBPos.h>
#include <protomol/output/OutputFinalXYZPos.h>
#include <protomol/output/OutputFinalXYZVel.h>
#include <protomol/output/OutputScreen.h>
#include <protomol/output/OutputXYZTrajectoryForce.h>
#include <protomol/output/OutputXYZTrajectoryPos.h>
#include <protomol/output/OutputXYZTrajectoryVel.h>
#include <protomol/output/OutputEnergies.h>
#include <protomol/output/OutputFAHGUI.h>
#include <protomol/output/OutputFAHFile.h>
#include <protomol/output/OutputScreen.h>
#include <protomol/output/OutputXTCTrajectory.h>
#include <protomol/output/OutputVelPosKE.h>
#include <protomol/output/OutputVelPosLS.h>
#include <protomol/output/OutputIonTemperatureLS.h>
#include <protomol/output/OutputIonSnapshot.h>

using namespace std;
using namespace ProtoMol;

void OutputModule::init(ProtoMolApp *app) {
  OutputFactory &f = app->outputFactory;

  f.registerExemplar(new OutputScreen());
  f.registerExemplar(new OutputDCDTrajectory());
  f.registerExemplar(new OutputDCDTrajectoryVel());
  f.registerExemplar(new OutputFinalPDBPos());
  f.registerExemplar(new OutputFinalXYZPos());
  f.registerExemplar(new OutputFinalXYZVel());
  f.registerExemplar(new OutputXYZTrajectoryForce());
  f.registerExemplar(new OutputXYZTrajectoryPos());
  f.registerExemplar(new OutputXYZTrajectoryVel());
  f.registerExemplar(new OutputEnergies());
#if defined (HAVE_GUI) || defined (HAVE_LIBFAH)
  f.registerExemplar(new OutputFAHGUI());
#endif
  f.registerExemplar(new OutputFAHFile());
  f.registerExemplar(new OutputXTCTrajectory());
  f.registerExemplar(new OutputVelPosLS()); 
  f.registerExemplar(new OutputVelPosKE()); 
  f.registerExemplar(new OutputIonTemperatureLS()); 
  f.registerExemplar(new OutputIonSnapshot()); 
}
