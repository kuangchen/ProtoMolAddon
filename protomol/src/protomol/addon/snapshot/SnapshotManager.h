#ifndef __SNAPSHOT_MANAGER_H
#define __SNAPSHOT_MANAGER_H

#include <protomol/addon/snapshot/SnapshotTimeQueue.h>
#include <fstream>
#include <vector>
#include <pair>
#include <string>

using std::vector;

namespace ProtoMolAddon {
  namespace Snapshot {
    
    template <class Storage> 
    class SnapshotManager {

    public:
      SnapshotManager(const ProtoMolApp *app, const string &fname) 
	: app(app) {

	ifstream is(fname);
	SnapshotTimeQueue tq;

	while (is >> tq) {
	  tq_list.push_back(tq);
	  storage_list.push_back(Storage(tq.size()));
	};
      }

      void Run(double t) {
	for (unsigned int i=0; i<tq_list.size(); i++)
	  if (tq_list[i].IsDue(t))
	    storage_list[i].Save(app);
      }

    private:
      const ProtoMolApp *app;
      vector<SnapshotTimeQueue> tq_list;
      vector<Storage> storage_list;
    };
  }
}
