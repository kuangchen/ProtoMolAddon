#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/GenericStorage.h>
#include <boost/format.hpp>
#include <iostream>

using std::cout;
using std::string;
using boost::format;
using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;


GenericStorage::GenericStorage(const string &fname) : fname(fname), current_frame(0) {}

void GenericStorage::SaveFrame(const ProtoMol::ProtoMolApp *app, double t) {
  current_frame++;
}

  
GenericStorage::~GenericStorage() {}

