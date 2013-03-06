#ifndef VESSELFORCES_MODULE_H
#define VESSELFORCES_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class VesselForcesModule : public Module {
  public:
    const std::string getName() const {return "VesselForces";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
    void registerForces(ProtoMolApp *app);
  };
};

#endif // VESSELFORCES_MODULE_H
