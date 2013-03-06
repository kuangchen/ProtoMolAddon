#ifndef NONBONDEDFULLFORCE_MODULE_H
#define NONBONDEDFULLFORCE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class NonbondedFullForceModule : public Module {
  public:
    NonbondedFullForceModule() {}
    const std::string getName() const {return "NonbondedFullForce";}
    void registerForces(ProtoMolApp *app);
  };
};

#endif // NONBONDEDFULLFORCE_MODULE_H
