#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/GenericStorage.h>
#include <boost/format.hpp>
#include <iostream>

using std::cout;
using std::string;
using boost::format;
using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;

GenericStorage::GenericStorage() {}

GenericStorage::GenericStorage(const string &fname, size_t total_frame_count) : 
  fname(fname), 
  total_frame_count(total_frame_count),
  save_time(0)
{}

void GenericStorage::Initialize(const ProtoMolApp *a) {
  app = a;
}

void GenericStorage::SaveFrame(double t) {
  save_time.push_back(t);
}

  
GenericStorage::~GenericStorage() {}

