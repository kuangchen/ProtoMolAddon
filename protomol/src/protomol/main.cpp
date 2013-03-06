#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>

#include <iostream>

using namespace std;
using namespace ProtoMol;

extern void moduleInitFunction(ModuleManager *);

int main(int argc, char *argv[]) {
  try {
    ModuleManager modManager;
    moduleInitFunction(&modManager);
    ProtoMolApp app(&modManager);

    if (!app.configure(argc, argv)) return 0;
    app.splash(cout);
    app.build();
    if ((int)app.config[InputDebug::keyword]) app.print(cout);

    while (app.step()) continue;
    app.finalize();

    return 0;

  } catch (const Exception &e) {
    cerr << "ERROR: " << e << endl;
  }

  return 1;
}

