#include "OutputModuleUnits.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/config/Configuration.h>
#include <protomol/factory/OutputFactory.h>

#include "../output/OutputScreenUnit.h"
#include "../output/OutputDCDTrajectoryPlus.h"

using namespace std;
using namespace ProtoMol;

void OutputModuleUnits::init(ProtoMolApp *app) {
  OutputFactory &f = app->outputFactory;

  //new scren for units
  f.unregisterExemplar("Screen");
  
  f.registerExemplar(new OutputScreenUnit());

  //new DCD for additional data hidden in comment
  f.unregisterExemplar("DCDFile");

  f.registerExemplar(new OutputDCDTrajectoryPlus());

}
