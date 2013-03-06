#ifndef HESSIANINTEGRATOR_MODULE_H
#define HESSIANINTEGRATOR_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class HessianIntegratorModule : public Module {
  public:
    const std::string getName() const {return "HessianIntegrator";}
    void init(ProtoMolApp *app);
  };
};

#endif // HESSIANINTEGRATOR_MODULE_H
