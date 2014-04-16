#include <protomol/addon/damping/DampingSpec.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

using namespace ProtoMolAddon::Damping;
using std::invalid_argument;
using std::cout;
using std::pair;
using std::make_pair;
using std::endl;
using std::cerr;

DampingSpec::DampingSpecEntry::DampingSpecEntry(double t_start, double t_end, double alpha):
  t_start(t_start), t_end(t_end), alpha(alpha) {
  CheckValidity();
}

void DampingSpec::DampingSpecEntry::CheckValidity() {
  if (t_start > t_end)
    throw invalid_argument("t_start > t_end");
  
  if (alpha < 0) 
    throw invalid_argument("alpha < 0");
}

namespace ProtoMolAddon { 
  namespace Damping {
 
    istream& operator>> (istream &is, DampingSpec::DampingSpecEntry &e) {
      try {
	is >> e.t_start >> e.t_end >> e.alpha;
	e.CheckValidity();
      }
      catch (const invalid_argument &e) {
	cerr << e.what() << endl;
      }

      return is;
    }

    ostream& operator<< (ostream &os, const DampingSpec::DampingSpecEntry &e) {
      os << "t_start = " << e.t_start
	 << "t_end = " << e.t_end
	 << "alpha = " << e.alpha << endl;

      return os;
    }
  }
}

DampingSpec::DampingSpec() : 
  entry_map() {}

istream& ProtoMolAddon::Damping::operator>> (istream &is, DampingSpec &spec) {
  string ion_name;
  DampingSpec::DampingSpecEntry entry;

  try {
    while (is >> ion_name >> entry) {
      pair<map<string, DampingSpec::DampingSpecEntry>::iterator, bool> ret;
      ret = spec.entry_map.insert(pair<string, DampingSpec::DampingSpecEntry>(ion_name, entry));

      if (!ret.second)
	cout << "Entry " << ion_name << " already exists, ignoring this entry " << entry << endl;
    }
  }

  catch (const invalid_argument &e) {
    cerr << e.what() << endl;
  }

  return is;
}

namespace ProtoMolAddon {
  namespace Damping {

    ostream& operator<< (ostream &os, const DampingSpec &spec) {
      for (const auto &kv : spec.entry_map) 
	os << kv.first << "\t" << kv.second << endl;

      return os;
    }
  }
}

void DampingSpec::Damp(ProtoMolApp *app, double now, double h) const  {

  for (unsigned int i=0; i<app->positions.size(); i++) {
    map<string, DampingSpec::DampingSpecEntry>::const_iterator it = entry_map.find(app->topology->atoms[i].name);

    if (it != entry_map.end()) {
      const DampingSpec::DampingSpecEntry &entry = it->second;

      if (now > entry.t_start && now < entry.t_end)  {
	app->velocities[i] = app->velocities[i] * exp(-entry.alpha * h);
	//app->positions[i] = app->positions[i] * exp(-entry.alpha * h);
      }
    }
  }
}
