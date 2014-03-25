#ifndef _TXTSTORAGE_H
#define _TXTSTORAGE_H

#include <memory>
#include <fstream>
#include <protomol/addon/snapshot/GenericStorage.h>

using std::shared_ptr;
using std::ofstream;

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {

    class TxtStorage : public GenericStorage {
    private:
      static string fname_pattern;
      static string separator;
      static size_t counter;

    public:
      TxtStorage();
      void SaveFrame(const ProtoMol::ProtoMolApp *app, double t);

    private:
      shared_ptr<ofstream> pf;
    };

  }
}


#endif
