#ifndef NONBONDEDINTERMITTENTFULLFORCE_MODULE_H
#define NONBONDEDINTERMITTENTFULLFORCE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class NonbondedIntermittentFullForceModule : public Module {
  public:
    const std::string getName() const {return "NonbondedIntermittentFullForce";}
    void registerForces(ProtoMolApp *app);
  };
};

#endif // NONBONDEDINTERMITTENTFULLFORCE_MODULE_H
