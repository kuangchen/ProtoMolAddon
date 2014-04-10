#ifndef __SNAPSHOT_MANAGER_H
#define __SNAPSHOT_MANAGER_H

#include <cassert>
#include <protomol/ProtoMolApp.h>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <string>

using ProtoMol::ProtoMolApp;

namespace ProtoMolAddon {
  namespace Snapshot {

    using namespace std;
    
    template <class TimeQueue, class Storage> 
    class SnapshotManager {
    private:
      typedef pair<TimeQueue, Storage> time_storage_pair;
      
    private:
      vector<time_storage_pair> time_storage;

    public:
      SnapshotManager() {}

      SnapshotManager(const string &fname) {
	
	ifstream is(fname);
	if (!is)
	  throw runtime_error("Cannot open file " + fname);

	TimeQueue tq;
	while (is >> tq) 
	  time_storage.push_back(time_storage_pair(tq, Storage(tq.Size())));
      }

      void Initialize(const ProtoMolApp *app) {
	for (pair<TimeQueue, Storage> &ts : time_storage) 
	  ts.second.Initialize(app);
      }
      
      void Run(double now) {
	for (pair<TimeQueue, Storage> &ts : time_storage) {
	  if (ts.first.IsDue(now)) {
	    double t = ts.first.PopFront();
	    ts.second.SaveFrame(t);
	  }
	}
      }
    
    };
  }
}
 
#endif 
