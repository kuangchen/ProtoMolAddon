#ifndef IO_MODULE_H
#define IO_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  declareInputValue(InputPositions, STRING, NOTEMPTY);
  declareInputValue(InputVelocities, STRING, NOTEMPTY);
  declareInputValue(InputPSF, STRING, NOTEMPTY);
  declareInputValue(InputPAR, STRING, NOTEMPTY);
  declareInputValue(InputPDBScaling, BOOL, NOCONSTRAINTS);
  declareInputValue(InputDihedralMultPSF, BOOL, NOCONSTRAINTS);
  declareInputValue(InputSCPISM, STRING, NOTEMPTY);

  declareInputValue(InputGromacsTopo, STRING, NOTEMPTY);
  declareInputValue(InputGromacsParamPath, STRING, NOTEMPTY);
  declareInputValue(InputGromacsGBParameterFile, STRING, NOTEMPTY);

  declareInputValue(InputGromacsTprFile, STRING, NOTEMPTY);

  class IOModule : public Module {
  public:
    const std::string getName() const {return "IO";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
    void read(ProtoMolApp *app);
  };
};

#endif // IO_MODULE_H
