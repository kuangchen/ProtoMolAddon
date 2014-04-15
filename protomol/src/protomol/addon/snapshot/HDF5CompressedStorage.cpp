#include <protomol/addon/snapshot/HDF5CompressedStorage.h>
#include <protomol/ProtoMolApp.h>
#include <iostream>
#include <algorithm>
#include <protomol/addon/Constants.h>
#include <H5DataSet.h>
#include <protomol/type/Vector3D.h>
#include <boost/format.hpp>
#include <cstring>

using std::string;
using boost::format;
using std::copy;
using std::endl;
using std::cerr;
using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::Snapshot;
using ProtoMol::Vector3D;

size_t HDF5CompressedStorage::counter(0);
string HDF5CompressedStorage::fname_pattern("snapshot_%d.hd5");

HDF5CompressedStorage::HDF5CompressedStorage() : 
  GenericStorage() {
  
}

HDF5CompressedStorage::HDF5CompressedStorage(size_t total_frame_count, 				     
					     unsigned int flags) 
try : GenericStorage((format(fname_pattern) % (counter++)).str(), total_frame_count), 
	file(fname.c_str(), flags)
	{}
catch (FileIException& e) {
  e.printError();
}

void HDF5CompressedStorage::Initialize(const ProtoMol::ProtoMolApp *app) {
  GenericStorage::Initialize(app);

  try {
    atom_count = app->positions.size();
    dataspace_dim[0] = total_frame_count;
    dataspace_dim[1] = atom_count;
    dataspace_dim[2] = 6;
    dataspace = DataSpace(3, dataspace_dim);    
  
    hsize_t chunk_dim[3]{total_frame_count/4, atom_count/2, 3};
    plist.setChunk(3, chunk_dim);
    plist.setDeflate(9);
    dataset = file.createDataSet("data", PredType::NATIVE_DOUBLE, dataspace, plist);
  }

  catch (FileIException &e) {
    e.printError();
  }

}

void HDF5CompressedStorage::Finalize() {
  try {
    // Record saved time, only if no frames has ever been recored
    if (!save_time.empty()) {
      hsize_t save_time_dim[1] {total_frame_count};
      DataSpace save_time_dataspace(1, save_time_dim);
      DataSet save_time_dataset = file.createDataSet("save time", 
						     PredType::NATIVE_DOUBLE, 
						     save_time_dataspace);

      save_time_dataset.write(save_time.data(), PredType::NATIVE_DOUBLE);
      save_time_dataset.close();
      save_time_dataspace.close();
 
      vector<ProtoMol::Atom>::const_iterator longest = max_element(app->topology->atoms.begin(),
								   app->topology->atoms.end(),
								   [](const ProtoMol::Atom &lhs, 
								      const ProtoMol::Atom &rhs) {
								     return lhs.name.size() < rhs.name.size();
								   });

      // Save atom name here

      vector<const char *> atom_name_list(atom_count);
      for (size_t i=0; i<atom_count; i++) 
	atom_name_list[i] = app->topology->atoms[i].name.c_str();

      hsize_t atom_name_dim[1] { atom_count };
      DataSpace atom_name_dataspace(1, atom_name_dim);
      StrType strdatatype(PredType::C_S1, H5T_VARIABLE);
      DataSet atom_name_dataset = file.createDataSet("atom name",
						     strdatatype,
						     atom_name_dataspace);


      atom_name_dataset.write(atom_name_list.data(), strdatatype);
      atom_name_dataset.close();
      atom_name_dataspace.close();

    }
    
  }
  catch (DataSetIException &e) {
    e.printError();
  }
  catch (DataSpaceIException &e) {
    e.printError();
  }

  catch (FileIException &e) {
    e.printError();
  }


}
  
HDF5CompressedStorage::~HDF5CompressedStorage() {
  try {
    file.close();

  } catch (FileIException &e) {
    e.printError();
  }
}

void HDF5CompressedStorage::SaveFrame(double t) {
  try {
    // Select the hyperspace 
    hsize_t start[3] {save_time.size(), 0, 0};
    hsize_t count[3] {1, 1, 1};
    hsize_t stride[3] {1, atom_count, 6};
    hsize_t block[3] {1, atom_count, 6};

    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
    
    hsize_t mspace_dim[1] {atom_count * 6};
    DataSpace mspace(1, mspace_dim);
    hsize_t startm[1]{0}, countm[1]{1}, stridem[1]{atom_count *6}, blockm[1]{atom_count*6};
    mspace.selectHyperslab(H5S_SELECT_SET, countm, startm, stridem, blockm );

    double *buffer = new double[atom_count * 6];
    Vector3D vel, pos;

    for (size_t i=0; i<atom_count; i++) {
      double *head = &(buffer[i * 6]);
      vel = app->velocities[i] * ToSI::velocity;
      pos = app->positions[i] * ToSI::position;
      copy(pos.c, pos.c+3, head);
      copy(vel.c, vel.c+3, head+3);
    }

    dataset.write(buffer, PredType::NATIVE_DOUBLE, mspace, dataspace);
    GenericStorage::SaveFrame(t);
  }
  

  catch( DataSpaceIException error ) {
    error.printError();
  }

  // catch failure caused by the H5File operations
  catch( AttributeIException error ) {
    error.printError();
  }


  // catch failure caused by the H5File operations
  catch( FileIException error ) {
    error.printError();
  }

  // catch failure caused by the DataSet operations
  catch( DataSetIException error ) {
    error.printError();
  }

}
