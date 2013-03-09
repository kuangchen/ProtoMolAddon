#include <string>
#include <sstream>
#include <iostream>

#include <protomol/addon/Snapshot.h>

using namespace ProtoMolAddon::IonSnapshot;
using namespace std;

Snapshot::Snapshot(const string& filename, const SnapshotConfig& config) :
  filename(filename),
  nextFrame(0),
  file(filename.c_str(), H5F_ACC_TRUNC),
  config(config) {
}

Snapshot::~Snapshot() {
  file.close();
}
  

void Snapshot::AddHeader() {
  hsize_t dim[] = {1};
  CompType header(sizeof(SnapshotConfig));
  header.insertMember("start", HOFFSET(SnapshotConfig, start), 
		      PredType::NATIVE_DOUBLE);
  header.insertMember("fps", HOFFSET(SnapshotConfig, fps), 
		      PredType::NATIVE_INT);
  header.insertMember("numFrame", HOFFSET(SnapshotConfig, numFrame), 
		      PredType::NATIVE_INT);
  header.insertMember("numAtom", HOFFSET(SnapshotConfig, numAtom), 
  		      PredType::NATIVE_INT);
  header.insertMember("numFrameperMM", HOFFSET(SnapshotConfig, numFrameperMM), 
		      PredType::NATIVE_INT);

 
  DataSpace dataspace(1, dim);
  DataSet dataset = file.createDataSet("/config", header, dataspace);
  
  dataset.write(&config, header);
  dataset.close();
  dataspace.close();

  return;
}

void Snapshot::AddFrame(Real time, Real *data) {
  if (!IsActive(time)) 
    return;


  ostringstream stream;
  stream << "/frame_" << nextFrame;
  hsize_t dims[2] = {config.numAtom, 14};
  DataSpace dataspace(2, dims);
  DataSet dataset = file.createDataSet(stream.str(), 
			       PredType::NATIVE_DOUBLE, 
			       dataspace);
  
  dataset.write(data, PredType::NATIVE_DOUBLE);  

  dataset.close();
  dataspace.close();
  
  nextFrame++;
  return;
}

bool Snapshot::IsActive(Real time) {
  //  return time>start && (nextFrame!=numFrame);

  return time > config.start && nextFrame != config.numFrame;
}
