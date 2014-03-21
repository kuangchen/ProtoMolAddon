#ifndef _TXTSTORAGE_H
#define _TXTSTORAGE_H

#include <iosfwd>
#include <protomol/addon/snapshot/GenericStorage.h>

using std::ofstream;

namespace ProtoMolAddon {
  namespace Snapshot {

    class TxtStorage : public GenericStorage {
      static string fname_pattern;

    public:
      TxtStorage(size_t n);
      ~TxtStorage();
      void Save(const ProtoMolApp *app);

    private:
      ofstream f;
    };

  }
}


#endif
