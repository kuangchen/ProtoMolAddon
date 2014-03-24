#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/GenericStorage.h>
#include <boost/format.hpp>

using std::string;
using boost::format;
using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;

unsigned int GenericStorage::counter(0);
string GenericStorage::fname_pattern("snapshot_%d.txt");

GenericStorage::GenericStorage() :
  id(counter++), fname((format(GenericStorage::fname_pattern) % id).str()) {}
  
GenericStorage::~GenericStorage() {}

