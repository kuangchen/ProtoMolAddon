#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>

#include <iomanip>
#include <iostream>

using namespace std;
using namespace ProtoMol;

extern void moduleInitFunction(ModuleManager *);

void splash(ostream &stream);

int main(int argc, char *argv[]) {
  try {
    ModuleManager modManager;
    moduleInitFunction(&modManager);
    ProtoMolApp app(&modManager);

    if (!app.configure(argc, argv)) return 0;
    splash(cout);
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

void splash(ostream &stream) {
  const int w = 16;
  const char* PROTOMOL_HR = "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+";

  stream
    << headerRow("ProtoCell") << endl
    << setw(w) << "Description: ";
  fillFormat(stream, "A rapid PROTOtyping CELLular dynamics object-oriented "
             "component based framework.", w, w);
  stream
#ifdef HAVE_PACKAGE_H
    << setw(w) << "Version: " << PACKAGE_VERSION << endl
    << setw(w) << "SVN revision: " << PACKAGE_REVISION << endl
    << setw(w) << "Repository: " << PACKAGE_SOURCE_REPO << endl
    << setw(w) << "Homepage: " << PACKAGE_HOMEPAGE << endl
    << setw(w) << "Report bugs to: " << PACKAGE_BUGREPORT << endl
    << setw(w) << "Compiler: " << PACKAGE_COMPILER << " "
    << PACKAGE_COMPILER_VERSION << endl
    << setw(w) << "Flags: " << PACKAGE_COMPILER_FLAGS << endl
    << setw(w) << "Extra libs: " << PACKAGE_COMPILER_LIBS << endl
    << setw(w) << "Built by: " << PACKAGE_BUILT_BY << endl
    << setw(w) << "Build platform: " << PACKAGE_PLATFORM << endl
#endif // HAVE_PACKAGE_H
    << setw(w) << "Build date: " <<  __DATE__ << ", " << __TIME__ << endl
#ifdef HAVE_LIBFAH
    << setw(w) << "Checksumming: "
    << "Enabled for Folding@Home file protection." << endl
#endif // HAVE_LIBFAH
#ifdef HAVE_PACKAGE_H
    << setw(w) << "Please cite: ";
  fillFormat(stream, PACKAGE_CITE, w, w);
  stream
#endif // HAVE_PACKAGE_H
    << PROTOMOL_HR << endl;
}

