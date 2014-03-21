#ifndef __SNAPSHOT_MANAGER_H
#define __SNAPSHOT_MANAGER_H

#include <protomol/ProtoMolApp.h>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <string>

using ProtoMol::ProtoMolApp;
using std::istream_iterator;
using std::copy;
using std::transform;
using std::vector;

namespace ProtoMolAddon {
  namespace Snapshot {
    
    template <class Storage> 
    class SnapshotManager {

    private:
      struct TimeQueueSpec {
	double t0;
	double dt;
	size_t n;
	size_t current;

	TimeQueueSpec(double t0=0, double dt=0, size_t n=0): 
	  t0(t0), dt(dt), n(n), current(0) {}

	void Pop() { current++; }

	friend istream& operator>> (istream &is, TimeQueueSpec &spec) {
	  is >> spec.t0 >> spec.dt >> spec.n;
	  spec.current = 0;
	  
	  return is;
	}
      };

    private:
      const ProtoMolApp *app;
      vector<TimeQueueSpec> tqspec_list;
      vector<Storage> storage_list;

    public:
      SnapshotManager(const ProtoMolApp *app = NULL) 
	: app(app) {}

      void Load(const string &fname) {
	tqspec_list.clear();
	storage_list.clear();

	ifstream is(fname);
	
	istream_iterator<TimeQueueSpec> begin(is), end;
	copy(begin, end, back_inserter(tqspec_list));

	storage_list.resize(tqspec_list.size());
      }

      void Run(double now) {
	for (size_t i=0; i<tqspec_list.size(); i++) {
	  TimeQueueSpec &spec(tqspec_list[i]);

	  if ( (spec.t0 + spec.current * spec.dt < now) && 
	       (now < spec.t0 + spec.n * spec.dt) ) {
	    storage_list[i].Save(app);
	    spec.current++;
	  }
	}

      }
    };
  }
}
 
#endif 
