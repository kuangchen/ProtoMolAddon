#ifndef COMMAND_LINE_MODULE_H
#define COMMAND_LINE_MODULE_H

#include <protomol/base/Module.h>

#include <string>
#include <vector>

namespace ProtoMol {
  class ProtoMolApp;

  class CommandLineModule : public Module {
    ProtoMolApp *app;

  public:
    CommandLineModule() : app(0) {}

    const std::string getName() const {return "CommandLine";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {
      return "Adds standard command line options.";
    }

    void getDependencies(module_deps_t &deps) const {
      deps.insert("Configuration");
    }

    void init(ProtoMolApp *app);

    int listTopologies(const std::vector<std::string> &args);
    int listIntegrators(const std::vector<std::string> &args);
    int listForces(const std::vector<std::string> &args);
    int listOutputs(const std::vector<std::string> &args);
  };
};

#endif // COMMAND_LINE_MODULE_H
