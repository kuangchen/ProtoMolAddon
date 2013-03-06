#ifndef TOPOLOGY_MODULE_H
#define TOPOLOGY_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <vector>
#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  declareInputValue(InputBoundaryConditions, STRING, NOTEMPTY);
  declareInputValue(InputCellManager, STRING, NOTEMPTY);

  class TopologyModule : public Module {
  public:
    const std::string getName() const {return "Topology";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {
      return "";
    }

    // Module interface
    void init(ProtoMolApp *app);
    void configure(ProtoMolApp *app);
  };
};

#endif // TOPOLOGY_MODULE_H
