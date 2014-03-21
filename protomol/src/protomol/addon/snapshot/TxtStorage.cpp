#include <fstream>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/snapshot/TxtStorage.h>

using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;
using namespace ProtoMolAddon::Constant;
using std::endl;

string TxtStorage::fname_pattern("snapshot_%d.txt");
string TxtStorage::separator("===================");

void TxtStorage::Save(const ProtoMol::ProtoMolApp *app) {
  size_t atom_count = app->positions.size();

  for (unsigned int i=0; i<atom_count; i++) 
    *pf << app->topology->atoms.at(i).name << " " 
	<< app->positions[i] * ToSI::position << " " 
	<< app->velocities[i] * ToSI::velocity << endl;
  
  *pf << TxtStorage::separator << endl;
}

