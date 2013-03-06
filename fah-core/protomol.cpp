#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/config.h>

#include <cbang/Exception.h>
#include <cbang/log/Logger.h>
#include <cbang/os/DynamicLibrary.h>

#include <fah/core/Core.h>
#include <fah/core/ExitCode.h>

#include <iostream>
#include <new> // For std::bad_alloc

#if HAVE_MKL
#include <omp.h>
#endif

using namespace std;
using namespace cb;
using namespace FAH;
using namespace ProtoMol;

extern void moduleInitFunction(ModuleManager *);

class ProtoMolCore : public Core {
public:
  ProtoMolCore() : Core("ProtoMol", CoreType::PROTOMOL) {}

  int init(int argc, char *argv[]) {
    options.add("steps-per-gen")->setType(Option::INTEGER_TYPE);

    int ret = Core::init(argc, argv);
    if (ret) return ret;

#if HAVE_MKL
    omp_set_num_threads(getNumProcs());
#endif

    // Add config file to args
    getArgs().push_front((char *)"protomol.conf");

    return 0;
  }


  int main(int argc, char *argv[]) {
    try {
#if defined(DEBUG) && !defined(_WIN32)
      ProtoMol::Debugger::initStackTrace(argv[0]);
#endif

      ModuleManager modManager;
      moduleInitFunction(&modManager);

      ProtoMolApp app(&modManager);

      app.splash(cout);

      // Load configuration options
      if (!app.load(argc, argv)) return 1;

      // Modify configuration
      // Add outputs FAHFile and FAHGUI
      app.config["FAHFile"] = "../current.xyz";
      app.config["FAHGUI"] = "ProtoMol";

      // Setup checkpointing
      app.config["Checkpoint"] = getCheckpointFile();
      app.config["CheckpointFreq"] = INT_MAX; // Disable

      // Configure and build
      app.configure();
      app.build();

      // Find CheckpointOutput
      OutputCheckpoint *oCheckpt = 0;
      OutputCollection::const_iterator it;
      const OutputCollection *outputs = app.outputs;
      for (it = outputs->begin(); it != outputs->end(); it++)
        if ((oCheckpt = dynamic_cast<OutputCheckpoint *>(*it))) break;

      if (!oCheckpt) THROW("Could not find OutputCheckpoint");

      // Set F@H info
      unsigned gen = getUnit().gen();
      int stepsPerGen = getOptions()["steps-per-gen"].toInteger();
      int firstStep = gen * stepsPerGen;
      int lastStep = firstStep + stepsPerGen;
      int frameSize = stepsPerGen < 100 ? 1 : stepsPerGen / 100;
      app.lastStep = lastStep; // Force correct last step
      setInfo(stepsPerGen, frameSize);

      // Check that the unit is still on track
      if (app.currentStep < firstStep || lastStep < app.currentStep) {
        LOG_ERROR("Invalid step " << app.currentStep << " not in ["
                  << firstStep << ", " << lastStep << "]");

        return ExitCode::BAD_CORE_FILES;
      }

      // Print configuration
      app.print(cout);

      do {
        // Update shared info file etc.
        step(app.currentStep - firstStep);
        
        if (shouldCheckpoint()) {
          oCheckpt->doIt(app.currentStep);
          checkpoint();
        }
      } while (!shouldQuit() && app.step(min(frameSize, 100)));
      
      oCheckpt->doIt(app.currentStep);
      app.finalize();
      checkpoint();

      return 0;

    } catch (const std::bad_alloc &e) {
      LOG_ERROR("ProtoMol std::bad_alloc");
      
    } catch (const ProtoMol::Exception &e) {
      LOG_ERROR("ProtoMol ERROR: " << e.getMessage());
    }

    getUnit().type() = CorePacketType::FAULTY;

    return ExitCode::UNKNOWN_ERROR;
  }
};


int main(int argc, char *argv[]) {
  DynamicLibrary::setEnabled(false); // Can crash in static builds
  return ProtoMolCore().run(argc, argv);
}
