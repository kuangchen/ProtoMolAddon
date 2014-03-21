#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/GenericStorage.h>
#include <boost/format.hpp>
#include <string>

using std::string;
using boost::format;
using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;

unsigned int GenericStorage::counter = 0;
string GenericStorage::fname_pattern("");

void GenericStorage::SetFilenamePattern(const string &pattern) {
  GenericStorage::fname_pattern = pattern;
}
  
GenericStorage::GenericStorage() :
  id(GenericStorage::counter++), 
  fname((format(GenericStorage::fname_pattern) % id).str()) {

}

