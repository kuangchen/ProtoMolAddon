#include <protomol/base/ModuleManager.h>
#include <protomol/base/Module.h>

#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/config/CommandLine.h>
#include <protomol/base/StringUtilities.h>

#include <iostream>

using namespace ProtoMol;
using namespace std;

void ModuleManager::add(Module *m) {
  if (find(m->getName()))
    THROW("Module '" + m->getName() + "' already registered");

  nameMap[m->getName()] = m;
  modules.insert(m);
}

void ModuleManager::remove(Module *m) {
  nameMap.erase(m->getName());
  modules.erase(m);
}

Module *ModuleManager::find(const string &name) {
  return nameMap[name].get();
}

void ModuleManager::init(ProtoMolApp *app) {
  modules_t::iterator it;

  // Check dependencies
  vector<string> missing;

  for (it = modules.begin(); it != modules.end(); it++) {
    Module::module_deps_t deps;
    (*it)->getDependencies(deps);

    Module::module_deps_t::iterator it2;
    for (it2 = deps.begin(); it2 != deps.end(); it2++)
      if (!nameMap[*it2]) missing.push_back((*it)->getName() + "->" + *it2);
  }

  if (missing.size()) {
    string errStr = "The following module dependencies were not satisfied:";

    for (unsigned int i = 0; i < missing.size(); i++)
      errStr += string(" ") + missing[i];

    THROW(errStr);
  }

  // Initialize modules
  for (it = modules.begin(); it != modules.end(); it++)
    (*it)->init(app);

  // Register command line option
  // TODO: This should probably go in a module.
  CommandLineOption::ActionBase *action =
    new CommandLineOption::Action<ModuleManager>(this,
                                                 &ModuleManager::listAction);

  app->cmdLine.add(0, "modules", action,
    "List all registered modules and exit.");
}

void ModuleManager::configure(ProtoMolApp *app) {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++)
    (*it)->configure(app);
}

void ModuleManager::read(ProtoMolApp *app) {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++) {
    (*it)->read(app);
  }
}

void ModuleManager::registerForces(ProtoMolApp *app) {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++)
    (*it)->registerForces(app);
}

void ModuleManager::postBuild(ProtoMolApp *app) {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++)
    (*it)->postBuild(app);
}

void ModuleManager::addModifiers(ProtoMolApp *app) {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++)
    (*it)->addModifiers(app);
}

int ModuleManager::listAction() {
  modules_t::iterator it;

  for (it = modules.begin(); it != modules.end(); it++) {
    const string name = (*it)->getName();
    const string help = (*it)->getHelp();

    cout << name;

    if (help != "") fillFormat(cout, help, name.length(), 40);
    else cout << endl;
  }

  return -1;
}

