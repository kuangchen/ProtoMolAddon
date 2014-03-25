#ifndef _HDF5_STORAGE_H
#define _HDF5_STORAGE_H

#include <protomol/addon/snapshot/GenericStorage.h>
#include <H5Cpp.h>
#include <H5File.h>

using namespace H5;

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {
    
    class HDF5Storage : public GenericStorage {
      static size_t counter;
      static string fname_pattern;


    public:
      HDF5Storage(unsigned int flags=H5F_ACC_TRUNC); 
      ~HDF5Storage(); 
      
      void SaveFrame(const ProtoMol::ProtoMolApp *app, double t);
      
    private:
      H5File file;
    };
  }
}

#endif
