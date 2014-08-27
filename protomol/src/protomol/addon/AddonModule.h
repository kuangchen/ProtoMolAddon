#ifndef _ADDON_MODULE_H
#define _ADDON_MODULE_H

#include <protomol/base/Module.h>
#include <string>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {

  class AddonModule : public ProtoMol::Module {
  public:
    const std::string getName() const {return "Addon";}
    const std::string getHelp() const {return "Addons for simulating AMO systems";}

    void init(ProtoMol::ProtoMolApp *app);
    void registerForces(ProtoMol::ProtoMolApp *app);
  };
}
  
#endif   
