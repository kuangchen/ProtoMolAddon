#ifndef _HDF5_COMPRESSED_STORAGE_H
#define _HDF5_COMPRESSED_STORAGE_H

#include <protomol/addon/snapshot/GenericStorage.h>
#include <H5Cpp.h>
#include <H5File.h>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {
    
    using namespace H5;

    class HDF5CompressedStorage : public GenericStorage {
      static size_t counter;
      static string fname_pattern;

    public:
      HDF5CompressedStorage();
      HDF5CompressedStorage(size_t total_frame_count, unsigned int flags=H5F_ACC_TRUNC); 

      ~HDF5CompressedStorage(); 

      void Initialize(const ProtoMol::ProtoMolApp *a);
      void Finalize();
      void SaveFrame(double t);
      
    private:
      size_t atom_count;
      H5File file;
      
      hsize_t dataspace_dim[3];
      DataSpace dataspace;
      DataSet dataset;
      DSetCreatPropList plist;
    };
  }
}

#endif
