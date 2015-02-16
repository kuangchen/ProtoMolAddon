#include <protomol/addon/snapshot/HDF5Storage.h>
#include <protomol/addon/snapshot/HDF5Storage.h>
#include <protomol/ProtoMolApp.h>
#include <H5DataSet.h>
#include <protomol/type/Vector3D.h>
#include <boost/format.hpp>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace ProtoMolAddon {
  namespace Snapshot {

    HDF5Storage::HDF5Storage() {}

    HDF5Storage::HDF5Storage(const HDF5Storage &other):
      const_ap_array_ptr(other.const_ap_array_ptr ?
			 new Util::ConstSIAtomProxyArray(*other.const_ap_array_ptr): NULL),
      file(other.file)
    {
    }
    
    HDF5Storage::HDF5Storage(const std::string &fname):
      file(fname, H5F_ACC_TRUNC)
    {
    }

    void HDF5Storage::Initialize(const ProtoMol::ProtoMolApp *app,
				 const std::queue<double> &unsaved_time) {
    
      const_ap_array_ptr.reset(new Util::ConstSIAtomProxyArray(app));
      size_t n_atom = const_ap_array_ptr->size();

      hsize_t trajectory_dataspace_dim[3];
      trajectory_dataspace_dim[0] = unsaved_time.size();
      trajectory_dataspace_dim[1] = n_atom;
      trajectory_dataspace_dim[2] = 6;

      trajectory_dataspace = DataSpace(3, trajectory_dataspace_dim);

      trajectory_dataset = file.createDataSet("trajectory",
					      PredType::NATIVE_DOUBLE,
					      trajectory_dataspace);
					      
    }

    void HDF5Storage::Save(const std::queue<double> &unsaved_time,
			   const std::vector<double> &saved_time) {
	
      size_t n_atom = const_ap_array_ptr->size();
      hsize_t start[3] {saved_time.size(), 0, 0};
      hsize_t count[3] {1, 1, 1};
      hsize_t stride[3] {1, n_atom, 6};
      hsize_t block[3] {1, n_atom, 6};

      trajectory_dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

      hsize_t mspace_dim[1] {n_atom * 6};
      DataSpace mspace(1, mspace_dim);
      hsize_t startm[1]{0};
      hsize_t countm[1]{1};
      hsize_t stridem[1]{n_atom *6};
      hsize_t blockm[1]{n_atom*6};
      
      mspace.selectHyperslab(H5S_SELECT_SET, countm, startm, stridem, blockm);

      double *buffer = new double[n_atom * 6];
      ProtoMol::Vector3D vel, pos;

      for (size_t i=0; i<n_atom; i++) {
	double *head = &(buffer[i * 6]);

	Util::ConstSIAtomProxy& const_ap = (*const_ap_array_ptr)[i];

	vel = const_ap.GetVelocity();
	pos = const_ap.GetPosition();

	std::copy(pos.c, pos.c+3, head);
	std::copy(vel.c, vel.c+3, head+3);
      }

      trajectory_dataset.write(buffer,
			       PredType::NATIVE_DOUBLE,
			       mspace,
			       trajectory_dataspace);
    }

    void HDF5Storage::Finalize(const std::queue<double> &unsaved_time,
			       const std::vector<double> &saved_time) {
      // Write save time
      hsize_t saved_time_dim[1] { saved_time.size() };
      DataSpace saved_time_dataspace(1, saved_time_dim);
      DataSet saved_time_dataset = file.createDataSet("time stamp", 
						      PredType::NATIVE_DOUBLE, 
						      saved_time_dataspace);

      saved_time_dataset.write(saved_time.data(), PredType::NATIVE_DOUBLE);
      saved_time_dataset.close();
      saved_time_dataspace.close();
      
      trajectory_dataset.close();
      trajectory_dataspace.close();

      // Write atom name
      size_t n_atom = const_ap_array_ptr->size();
      std::vector<const char *> atom_name_list(n_atom);
      for (size_t i=0; i<n_atom; i++) 
	atom_name_list[i] = (*const_ap_array_ptr)[i].GetName().c_str();

      hsize_t atom_name_dim[1] { n_atom };
      DataSpace atom_name_dataspace(1, atom_name_dim);
      StrType strdatatype(PredType::C_S1, H5T_VARIABLE);
      DataSet atom_name_dataset = file.createDataSet("atom name",
						     strdatatype,
						     atom_name_dataspace);

      atom_name_dataset.write(atom_name_list.data(), strdatatype);
      atom_name_dataset.close();
      atom_name_dataspace.close();

      file.close();
    }

    
  }
}
