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
      vector<TimeQueue> time_queue;
      vector<Storage> storage;

    public:
      SnapshotManager() : time_queue(0), storage(0) {}

      SnapshotManager(const string &fname) {
	ifstream is(fname);
	if (!is)
	  throw runtime_error("Cannot open file " + fname);

	copy(istream_iterator<TimeQueue>(is), istream_iterator<TimeQueue>(), back_inserter(time_queue));

	for (int i=0; i<time_queue.size(); i++)
	  storage.push_back(Storage(time_queue[i].Size()));
      }

      void Initialize(const ProtoMolApp *app) {
	for (auto &s: storage)
	  s.Initialize(app);
      }
      
      void Run(double now) {
	for (int i=0; i<time_queue.size(); i++) 
	  if (time_queue[i].IsDue(now)) {
	    double t = time_queue[i].PopFront();
	    storage[i].SaveFrame(t);
	  }
      }
    
    };
  }
}
 
#endif 
