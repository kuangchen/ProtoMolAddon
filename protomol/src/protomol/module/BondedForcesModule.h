#ifndef BONDEDFORCES_MODULE_H
#define BONDEDFORCES_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class BondedForcesModule : public Module {
  public:
    const std::string getName() const {return "BondedForces";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
    void registerForces(ProtoMolApp *app);
  };
};

#endif // BONDEDFORCES_MODULE_H
