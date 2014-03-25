#include <fstream>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/TxtStorage.h>
#include <boost/format.hpp>

using ProtoMol::ProtoMolApp;
using boost::format;
using namespace ProtoMolAddon::Snapshot;
using namespace ProtoMolAddon::Constant;
using std::endl;

size_t TxtStorage::counter(0);
string TxtStorage::fname_pattern("snapshot_%d.txt");
string TxtStorage::separator("===================");

TxtStorage::TxtStorage()
  : GenericStorage((boost::format(fname_pattern) % (counter++)).str()), 
    pf(new ofstream(fname)) {

}

void TxtStorage::SaveFrame(const ProtoMol::ProtoMolApp *app, double t) {
  size_t atom_count = app->positions.size();
  
  *pf << t << endl;

  for (unsigned int i=0; i<atom_count; i++) 
    *pf << app->topology->atoms.at(i).name << " " 
	<< app->positions[i] * ToSI::position << " " 
	<< app->velocities[i] * ToSI::velocity << endl;
  
  *pf << TxtStorage::separator << endl;
}
