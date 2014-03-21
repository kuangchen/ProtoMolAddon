#include <protomol/addon/snapshot/TxtStorage.h>
#include <fstream>
#include <iostream>
#include <string>
#include <protomol/addon/Constants.h>
#include <protomol/ProtoMolApp.h>

using ProtoMol::ProtoMolApp;
using namespace ProtoMolAddon::Snapshot;
using namespace ProtoMolAddon::Constant;
using std::endl;

string ProtoMolAddon::Snapshot::TxtStorage::fname_pattern("snapshot_%d.txt");

TxtStorage::TxtStorage(size_t n):
  GenericStorage(),
  f(fname) {
}

void TxtStorage::Save(const ProtoMolApp *app) {
  size_t atom_count = app->positions.size();

  for (unsigned int i=0; i<atom_count; i++) {
    f << app->topology->atoms.at(i).name << " " 
      << app->positions[i] * ToSI::position << " " 
      << app->velocities[i] * ToSI::velocity << endl;
  }

  f << "============================================" << endl;
}
