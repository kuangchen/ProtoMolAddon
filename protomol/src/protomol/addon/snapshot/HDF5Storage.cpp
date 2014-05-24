#include <protomol/addon/snapshot/HDF5Storage.h>
#include <protomol/ProtoMolApp.h>
#include <iostream>
#include <algorithm>
#include <protomol/addon/util/SIAtomProxy.h>
#include <H5DataSet.h>
#include <protomol/type/Vector3D.h>
#include <boost/format.hpp>
#include <cstring>

using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Snapshot;
using ProtoMol::Vector3D;

size_t HDF5Storage::file_name_counter(0);

string HDF5Storage::file_name_pattern("snapshot_%d.hd5");

void HDF5Storage::SetFileNamePattern(const std::string &pattern) { file_name_pattern = pattern; }

HDF5Storage::HDF5Storage() {}

HDF5Storage::HDF5Storage(const TimeQueue &tq, unsigned int flags) :
  tq(tq),
  file((boost::format(file_name_pattern) % (file_name_counter++)).str(), flags)
{}


void HDF5Storage::Initialize(const ProtoMol::ProtoMolApp *app_) {
  app = app_;

  atom_count = app->positions.size();
  dataspace_dim[0] = tq.Size();
  dataspace_dim[1] = atom_count;
  dataspace_dim[2] = 6;
  dataspace = H5::DataSpace(3, dataspace_dim);

  size_t d = atom_count == 1 ? atom_count : atom_count / 2;
  
  hsize_t chunk_dim[3]{tq.Size()/4, d, 3};
  plist.setChunk(3, chunk_dim);
  plist.setDeflate(9);
  dataset = file.createDataSet("data", PredType::NATIVE_DOUBLE, dataspace, plist);
}

void HDF5Storage::Finalize() {
  // Record saved time, only if no frames has ever been recored
  if (tq.SavedTime().empty()) return;

  hsize_t save_time_dim[1] { tq.SavedTime().size() };
  DataSpace save_time_dataspace(1, save_time_dim);
  DataSet save_time_dataset = file.createDataSet("save time", 
						 PredType::NATIVE_DOUBLE, 
						 save_time_dataspace);

  std::vector<double> saved_time(tq.SavedTime().begin(), tq.SavedTime().end());
  save_time_dataset.write(saved_time.data(), PredType::NATIVE_DOUBLE);
  save_time_dataset.close();
  save_time_dataspace.close();
 
  // Save atom name here

  std::vector<const char *> atom_name_list(atom_count);
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
  

HDF5Storage::~HDF5Storage() {
  file.close();
}

void HDF5Storage::SaveFrame(double now) {
  if (!tq.IsDue(now)) return;

    // Select the hyperspace 
  hsize_t start[3] {tq.SavedTime().size(), 0, 0};
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

    Util::ConstSIAtomProxy atom(app, i);
      
    vel = atom.GetVelocity();
    pos = atom.GetPosition();

    std::copy(pos.c, pos.c+3, head);
    std::copy(vel.c, vel.c+3, head+3);
  }

  dataset.write(buffer, PredType::NATIVE_DOUBLE, mspace, dataspace);
  delete buffer;

  tq.Pop();
}
