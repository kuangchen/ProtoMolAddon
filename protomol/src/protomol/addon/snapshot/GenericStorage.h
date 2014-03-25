#ifndef _GENERIC_STORAGE_H
#define _GENERIC_STORAGE_H

#include <string>

namespace ProtoMol {
  class ProtoMolApp;
}

using std::string;

namespace ProtoMolAddon {
  namespace Snapshot {

    class GenericStorage {
    protected:
      string fname;
      size_t current_frame;

    public:
      GenericStorage(const string &fname);
	
      virtual ~GenericStorage() = 0; 
      virtual void SaveFrame(const ProtoMol::ProtoMolApp *app, double t) = 0;

    };

  }
}

#endif
