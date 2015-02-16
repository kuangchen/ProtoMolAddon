#ifndef __HDF5_STORAGE_H
#define __HDF5_STORAGE_H

#include <H5Cpp.h>
#include <H5File.h>

#include <memory>
#include <string>
#include <iostream>
#include <queue>
#include <protomol/addon/util/ConstSIAtomProxy.h>
#include <protomol/addon/util/ConstSIAtomProxyArray.h>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxyArray;
    class ConstSIAtomProxy;
  }
  
  namespace Snapshot {
    using namespace H5;

    class HDF5Storage {
    private:
      std::unique_ptr<Util::ConstSIAtomProxyArray> const_ap_array_ptr;
      H5File file;

      DataSpace trajectory_dataspace;
      DataSet trajectory_dataset;
      
    public:
      HDF5Storage();
      HDF5Storage(const HDF5Storage &other);
      HDF5Storage(const std::string &fname);
      
      void Initialize(const ProtoMol::ProtoMolApp *app,
		      const std::queue<double> &unsaved_time);

      void Save(const std::queue<double> &unsaved_time,
		const std::vector<double> &saved_time);

      void Finalize(const std::queue<double> &unsaved_time,
		    const std::vector<double> &saved_time);

      static std::string GetName() { return "HDF5"; }
      static std::string GetParameterName() { return "-HDF5-spec"; }
    };
  }
}

#endif 
