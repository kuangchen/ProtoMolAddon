#include <protomol/addon/snapshot/HDF5Storage.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/addon/util/SIAtomProxyArray.h>
#include <H5DataSet.h>
#include <protomol/type/Vector3D.h>
#include <boost/format.hpp>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace ProtoMolAddon {
  namespace Snapshot {

    size_t HDF5Storage::file_name_counter(0);

    std::string HDF5Storage::file_name_pattern("snapshot_%d.hd5");

    void HDF5Storage::SetFileNamePattern(const std::string &pattern) { file_name_pattern = pattern; }

    HDF5Storage::HDF5Storage() {}

    HDF5Storage::HDF5Storage(const TimeQueue &tq, unsigned int flags) :
      tq(tq),
      file((boost::format(file_name_pattern) % (file_name_counter++)).str(), flags)
    {}


    void HDF5Storage::Initialize(const ProtoMolApp *app) {
      ap_array_ptr.reset(new Util::ConstSIAtomProxyArray(app));

      n_atom = ap_array_ptr->size();

      tr_dim[0] = tq.size();
      tr_dim[1] = n_atom;
      tr_dim[2] = 6;
      tr_dspace = DataSpace(3, tr_dim);

      // Use compression
      size_t d = n_atom == 1 ? n_atom : n_atom / 2;
      hsize_t chunk_dim[3]{tq.size()/4, d, 3};
      tr_plist.setChunk(3, chunk_dim);
      tr_plist.setDeflate(9);
      
      tr_dset = file.createDataSet("data",
				   PredType::NATIVE_DOUBLE,
				   tr_dspace,
				   tr_plist);

      // Write atom's name here
      std::vector<const char *> atom_name_list;
      for (auto &ap: *ap_array_ptr)
	atom_name_list.push_back(ap.GetName().c_str());

      hsize_t name_dim[1] { n_atom };
      DataSpace name_dspace(1, name_dim);
      StrType strdatatype(PredType::C_S1, H5T_VARIABLE);
      DataSet name_dset = file.createDataSet("atom name",
					     strdatatype,
					     name_dspace);

      name_dset.write(atom_name_list.data(), strdatatype);
      name_dset.close();
      name_dspace.close();
    }

    HDF5Storage::~HDF5Storage() { file.close(); }
    
    void HDF5Storage::Finalize() {
      // Record saved time, only if no frames has ever been recored
      size_t saved_time_dim = tq.get_init_len() - tq.size();
      
      if (saved_time_dim==0) return;

      hsize_t save_time_dim[1] { saved_time_dim };
      DataSpace save_time_dspace(1, saved_time_dim);
      DataSet save_time_dset = file.createDataSet("save time", 
						  PredType::NATIVE_DOUBLE, 
						  save_time_dspace);

      std::vector<double> saved_time(tq.cbegin(), tq.current());
      save_time_dset.write(saved_time.data(), PredType::NATIVE_DOUBLE);
      save_time_dset.close();
      save_time_dspace.close();
    }
  
    void HDF5Storage::SaveFrame(double now) {
      if (tq.front() < now) return;
      
      // Select the hyperspace
      // tq.saved_size(), 0, 0 to tq.saved_size()+1, n_atom, 6
      hsize_t offset[3] {tq.saved_size(), 0, 0};
      hsize_t count[3] {1, n_atom, 6};
      tr_dspace.selectHyperslab(H5S_SELECT_SET, count, offset, NULL, NULL);
    
      hsize_t buffer_dim[1] {n_atom * 6};
      DataSpace buffer_dspace(1, buffer_dim);
      hsize_t buffer_offset[1]{0};
      hsize_t buffer_count[1]{ n_atom * 6};
      buffer_dspace.selectHyperslab(H5S_SELECT_SET, buffer_count, buffer_offset);

      double *buffer = new double[n_atom * 6];

      for (auto &ap: *ap_array_ptr) {
	double *head = &(buffer[ap.GetID()*6]);
	const Vector3D &pos = ap.GetPosition();
	const Vector3D &vel = ap.GetVelocity();

	std::copy(pos.c, pos.c+3, head);
	std::copy(vel.c, vel.c+3, head+3);
      }

      tr_dset.write(buffer, PredType::NATIVE_DOUBLE, buffer_dspace, tr_dspace);
      delete buffer;

      tq.pop_front();
      tr_dspace.selectNone();
    }

  }
}
