#ifndef NONBONDEDCUTOFFFORCE_MODULE_H
#define NONBONDEDCUTOFFFORCE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class NonbondedCutoffForceModule : public Module {
  public:
    const std::string getName() const {return "NonbondedCutoffForce";}
    void registerForces(ProtoMolApp *app);
  };
};

#endif // NONBONDEDCUTOFFFORCE_MODULE_H
