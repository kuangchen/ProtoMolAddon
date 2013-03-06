#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/config/InputValue.h>
#include <protomol/base/Report.h>

#include <iostream>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

extern void moduleInitFunction(ModuleManager *);

//define our config file input values here
namespace ProtoMol{
  declareInputValue(InputDcd, STRING, NOTEMPTY);
  defineInputValue(InputDcd, "inputDcdFile");
}

int main(int argc, char *argv[]) {  

  if(argc < 2){
    report << plain << "Require .conf file name." << endr;
    return 0;
  }

  // dcd file
  DCDTrajectoryReader in;
  try {
    TimerStatistic::timer[TimerStatistic::WALL].start();
    ModuleManager modManager;
    moduleInitFunction(&modManager);
    ProtoMolApp app(&modManager);

    //register our config values
    InputDcd::registerConfiguration(&(app.config));

    //pass config file name
    vector<string> args;
    args.push_back("ProtoMol");
    args.push_back(argv[1]);
    if (!app.configure(args)) return 0;

    app.splash(cout);
    app.build();
    if ((int)app.config[InputDebug::keyword]) app.print(cout);

    //open dcd file
    if(in.open(app.config[InputDcd::keyword])){
      try{
        in >> app.positions;
      }catch(const Exception &e){ //end if read error
        app.currentStep = app.lastStep;
      }
    }else{    
      report << plain << "Input file '" << app.config[InputDcd::keyword] << "' not defined." << endr;
      return 0;
    }
    
    //loop until end of dcd OR number of steps
    while(app.currentStep < app.lastStep){
      app.outputs->run(app.currentStep);

      int inc = app.outputs->getNext() - app.currentStep;
      inc = std::min(app.lastStep, app.currentStep + inc) - app.currentStep;
      app.currentStep += inc;

      try{
        in >> app.positions;
      }catch(const Exception &e){ //end if read error
        app.currentStep = app.lastStep;
      }

    }
    //
    app.finalize();
    TimerStatistic::timer[TimerStatistic::WALL].stop();
    return 0;

  } catch (const Exception &e) {
    cerr << "ERROR: " << e << endl;
  }

  return 1;
}



