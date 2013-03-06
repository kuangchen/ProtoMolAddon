#ifndef TEMPLATE_MODULE_H
#define TEMPLATE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class TemplateModule : public Module {
  public:
    const std::string getName() const {return "Template";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
  };
};

#endif // TEMPLATE_MODULE_H
