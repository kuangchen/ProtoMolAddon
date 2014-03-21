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

    public:
      TxtStorage(): GenericStorage(), pf(new ofstream(fname)) {}

      void Save(const ProtoMol::ProtoMolApp *app);

    private:
      shared_ptr<ofstream> pf;
    };

  }
}


#endif
