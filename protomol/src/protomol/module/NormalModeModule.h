#ifndef NORMALMODE_MODULE_H
#define NORMALMODE_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  declareInputValue(InputEigenVectors, STRING, NOTEMPTY);
  declareInputValue(InputEigTextFile, STRING, NOTEMPTY);
  declareInputValue(InputEigenValues, STRING, NOTEMPTY);

  class NormalModeModule : public Module {

  public:
    NormalModeModule() {}

    const std::string getName() const {return "NormalMode";}
    int getPriority() const {return 5;} // Must be after IOModule

    void init(ProtoMolApp *app);
    void read(ProtoMolApp *app);
  };
};

#endif // NORMALMODE_MODULE_H
