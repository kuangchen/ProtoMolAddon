#ifndef INTEGRATOR_OPEMM_MODULE_H
#define INTEGRATOR_OPEMM_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class IntegratorOpenMMModule : public Module {
  public:
    const std::string getName() const {return "IntegratorOpenMM";}
    void init(ProtoMolApp *app);
  };
};

#endif // INTEGRATOR_OPEMM_MODULE_H
