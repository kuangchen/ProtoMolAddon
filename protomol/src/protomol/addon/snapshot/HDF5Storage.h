#ifndef _HDF5_STORAGE_H
#define _HDF5_STORAGE_H

#include <string>
#include <H5Cpp.h>
#include <H5File.h>
#include <protomol/addon/snapshot/TimeQueue.h>
#include <protomol/addon/util/ConstSIAtomProxyArray.h>
#include <memory>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {
    
    using namespace H5;
    using namespace ProtoMol;
    
    class HDF5Storage {
    private:
      static size_t file_name_counter;
      static std::string file_name_pattern;
      
    public:
      static void SetFileNamePattern(const std::string &pattern);

      static std::string GetName() { return "HDF5Storage"; }
      static std::string GetParameterName() { return "-hd5-storage-spec"; }

    private:
      TimeQueue tq;
      std::unique_ptr<Util::ConstSIAtomProxyArray> ap_array_ptr;
      
      H5File file;
      size_t n_atom;

      hsize_t tr_dim[3];
      DataSpace tr_dspace;
      DataSet tr_dset;
      DSetCreatPropList tr_plist;
      
      //      hsize_t dataspace_dim[3];
      //DataSpace dataspace;
      //DataSet dataset;
      //DSetCreatPropList plist;

    public:
      
      HDF5Storage();
      ~HDF5Storage();
      HDF5Storage(const TimeQueue &tq, unsigned int flags=H5F_ACC_TRUNC); 
      HDF5Storage(HDF5Storage &&) = default;
      void Initialize(const ProtoMolApp *app);
      void Finalize();
      void SaveFrame(double t);
    };
  }
}

#endif
