#ifndef MAIN_MODULE_H
#define MAIN_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

namespace ProtoMol {
  class ProtoMolApp;
  
  declareInputValue(InputNumsteps, INT, NOTNEGATIVE);
  declareInputValue(InputFirststep, INT, NOCONSTRAINTS);
  declareInputValue(InputRealfirststep, INT, NOCONSTRAINTS);
  declareInputValue(InputSeed, INT, NOTNEGATIVE);
  declareInputValue(InputRandomType, INT, NOTNEGATIVE);
  declareInputValue(InputDebug, INT, NOCONSTRAINTS);
  declareInputValue(InputIntegrator, INTEGRATOR, NOTEMPTY);
  declareInputValue(InputReducedImage, BOOL, NOCONSTRAINTS);
  declareInputValue(InputTemperature, REAL, NOTNEGATIVE);
  declareInputValue(InputDoSCPISM, BOOL, NOCONSTRAINTS);
  declareInputValue(InputMolVirialCalc, BOOL, NOCONSTRAINTS);
  declareInputValue(InputVirialCalc, BOOL, NOCONSTRAINTS);
  declareInputValue(InputOutputfreq, INT, NOTNEGATIVE);
  declareInputValue(InputOutput, BOOL, NOCONSTRAINTS);
  declareInputValue(InputMinimalImage, BOOL, NOCONSTRAINTS);
  declareInputValue(InputDoGBSAObc, BOOL, NOCONSTRAINTS);
  declareInputValue(InputDebugLimit, INT, NOCONSTRAINTS);


  class MainModule : public Module {
  public:
    const std::string getName() const {return "Main";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {
      return "System default keywords and setup.";
    }

    void init(ProtoMolApp *app);
    void configure(ProtoMolApp *app);
    void postBuild(ProtoMolApp *app);
  };
};

#endif // MAIN_MODULE_H
