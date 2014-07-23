#ifndef _HDF5_STORAGE_H
#define _HDF5_STORAGE_H

#include <string>
#include <H5Cpp.h>
#include <H5File.h>
#include <protomol/addon/snapshot/TimeQueue.h>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Snapshot {
    
    using namespace H5;

    class HDF5Storage {
    private:
      static size_t file_name_counter;
      static std::string file_name_pattern;
      
    public:
      static void SetFileNamePattern(const std::string &pattern);
      static std::string GetName() { return "HDF5Storage"; }

    private:
      TimeQueue tq;
      const ProtoMol::ProtoMolApp *app;

      H5File file;

      size_t atom_count;      
      hsize_t dataspace_dim[3];
      DataSpace dataspace;
      DataSet dataset;
      DSetCreatPropList plist;

    public:
      HDF5Storage();
      HDF5Storage(const TimeQueue &tq, unsigned int flags=H5F_ACC_TRUNC); 

      ~HDF5Storage(); 

      void Initialize(const ProtoMol::ProtoMolApp *a);
      void Finalize();
      void SaveFrame(double t);

    };
  }
}

#endif
