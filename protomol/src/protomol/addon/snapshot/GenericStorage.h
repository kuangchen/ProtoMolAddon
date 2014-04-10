#ifndef _GENERIC_STORAGE_H
#define _GENERIC_STORAGE_H

#include <string>

namespace ProtoMol {
  class ProtoMolApp;
}


namespace ProtoMolAddon {
  namespace Snapshot {

    using std::string;

    class GenericStorage {
    protected:
      string fname;
      size_t current_frame;
      size_t total_frame_count;
      const ProtoMol::ProtoMolApp *app;

    public:
      
      GenericStorage(const string &fname, size_t total_frame_count);
	
      virtual ~GenericStorage() = 0; 
      void Initialize(const ProtoMol::ProtoMolApp *a);
      virtual void SaveFrame(double t) = 0;

    };
  }
}

#endif
