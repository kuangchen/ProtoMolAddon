#ifndef NONBONDEDFULLELECTROSTATICFORCE_MODULE_H
#define NONBONDEDFULLELECTROSTATICFORCE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class NonbondedFullElectrostaticForceModule : public Module {
  public:
    NonbondedFullElectrostaticForceModule() {}
    const std::string getName() const
    {return "NonbondedFullElectrostaticForce";}
    void registerForces(ProtoMolApp *app);
  };
};

#endif // NONBONDEDFULLELECTROSTATICFORCE_MODULE_H
