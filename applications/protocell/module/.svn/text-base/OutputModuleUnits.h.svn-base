#ifndef OUTPUT_MODULE_UNITS_H
#define OUTPUT_MODULE_UNITS_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;
  
  class OutputModuleUnits : public Module {
  public:
    const std::string getName() const {return "OutputUnits";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
  };
};

#endif // OUTPUT_MODULE_UNITS_H
