#ifndef ADDON_MODULE_H
#define ADDON_MODULE_H

#include <protomol/base/Module.h>
#include <string>

namespace ProtoMolAddon {

  class AddonModule : public ProtoMol::Module {
  public:
    const std::string getName() const {return "Addon";}
    void init(ProtoMol::ProtoMolApp *app);
    const std::string getHelp() const {return "Addons for simulating AMO systems";}
  };

}
  
#endif   
