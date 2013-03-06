#ifndef MODULE_MANAGER_H
#define MODULE_MANAGER_H

#include <string>
#include <set>
#include <map>

#include <protomol/base/Module.h>
#include <protomol/base/SmartPointer.h>

namespace ProtoMol {
  class ProtoMolApp;

  class ModuleManager {
    typedef std::map<std::string, SmartPointer<Module> > nameMap_t;
    nameMap_t nameMap;

    typedef std::set<Module *, moduleLess> modules_t;
    modules_t modules;

  public:
    void add(Module *m);
    void remove(Module *m);
    Module *find(const std::string &name);

    void init(ProtoMolApp *app);
    void configure(ProtoMolApp *app);
    void read(ProtoMolApp *app);
    void registerForces(ProtoMolApp *app);
    void postBuild(ProtoMolApp *app);
    void addModifiers(ProtoMolApp *app);

    int listAction();
  };
};

#endif // MODULE_MANAGER_H
