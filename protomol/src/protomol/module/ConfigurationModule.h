#ifndef CONFIGURATION_MODULE_H
#define CONFIGURATION_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <vector>
#include <string>

namespace ProtoMol {
  class ProtoMolApp;
  class Configuration;

  declareInputValue(InputConfig, STRING, NOTEMPTY);

  class ConfigurationModule : public Module {
    Configuration *config;

  public:
    const std::string getName() const {return "Configuration";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {
      return "Configuration system.";
    }

    // From Module
    void init(ProtoMolApp *app);
    using Module::configure; // To avoid compiler warning/error

    int listKeywords();
    int configure(const std::vector<std::string> &args);
  };
};

#endif // CONFIGURATION_MODULE_H
