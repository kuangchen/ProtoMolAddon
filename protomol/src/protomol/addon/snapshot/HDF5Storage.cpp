#include <protomol/addon/snapshot/HDF5Storage.h>
#include <protomol/ProtoMolApp.h>
#include <iostream>
#include <algorithm>
#include <protomol/addon/Constants.h>
#include <H5DataSet.h>
#include <protomol/type/Vector3D.h>
#include <boost/format.hpp>


using std::string;
using boost::format;
using std::copy;
using std::endl;
using std::cerr;
using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::Snapshot;
using ProtoMol::Vector3D;

size_t HDF5Storage::counter(0);
string HDF5Storage::fname_pattern("snapshot_%d.hd5");


HDF5Storage::HDF5Storage(unsigned int flags)
try : GenericStorage((format(fname_pattern) % (counter++)).str()), 
	file(fname.c_str(), flags) {
}

catch (FileIException& e) {
  e.printError();
}
  
HDF5Storage::~HDF5Storage() {
  try {
    file.close();
  } catch (FileIException &e) {
    e.printError();
  }
}

void HDF5Storage::SaveFrame(const ProtoMol::ProtoMolApp *app, double t) {

  try {
  
    size_t atom_count(app->positions.size());
    hsize_t dataset_dim[2] = { atom_count, 6 };
  
    // // Write the velocity and position
    DataSpace dataspace(2, dataset_dim);
  
    DataSet dataset = file.createDataSet((boost::format("frame_%d")% current_frame).str().c_str(), 
     					 PredType::NATIVE_DOUBLE, dataspace);
    
    double *data = new double[atom_count * 6];
    Vector3D vel, pos;

    for (size_t i=0; i<atom_count; i++) {
      double *head = &(data[i * 6]);
      vel = app->velocities[i] * ToSI::velocity;
      pos = app->positions[i] * ToSI::position;

      copy(pos.c, pos.c+3, head);
      copy(vel.c, vel.c+3, head+3);
    }

    dataset.write(data, PredType::NATIVE_DOUBLE);

    // Write the attribute
    double attr_data[1] = { t };
    hsize_t attr_dim[1] = { 1 };
    DataSpace attr_dataspace(1, attr_dim);
    Attribute attribute = dataset.createAttribute("Time", 
    						  PredType::NATIVE_DOUBLE,
    						  attr_dataspace,
    						  PropList::DEFAULT );

    attribute.write(PredType::NATIVE_DOUBLE, attr_data);

    dataset.close();
    attribute.close();
    delete data;
    GenericStorage::SaveFrame(app, t);
  }

// catch failure caused by the H5File operations
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
