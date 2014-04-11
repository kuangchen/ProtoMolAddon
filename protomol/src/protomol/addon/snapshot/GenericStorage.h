#ifndef _GENERIC_STORAGE_H
#define _GENERIC_STORAGE_H

#include <string>
#include <vector>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {

    using std::string;
    using std::vector;

    class GenericStorage {
    protected:
      string fname;
      size_t total_frame_count;

      const ProtoMol::ProtoMolApp *app;
      vector<double> save_time;

    public:
      GenericStorage();
      GenericStorage(const string &fname, size_t total_frame_count);
	
      virtual ~GenericStorage() = 0; 

      void Initialize(const ProtoMol::ProtoMolApp *a);
      void Finalize();

      virtual void SaveFrame(double t) = 0;
    };
  }
}

#endif
