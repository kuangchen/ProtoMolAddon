#ifndef __SNAPSHOT_TIMEQUEUE_H
#define __SNAPSHOT_TIMEQUEUE_H

#include <iosfwd>
#include <queue>

using std::istream;
using std::queue;

namespace ProtoMolAddon {
  namespace Snapshot {

    class SnapshotTimeQueue : private queue<double> {
    public:
      SnapshotTimeQueue();

    private:
      void Reset(double t0, double dt, size_t n);
      inline bool IsDue(double t);
      
      friend istream& operator>> (istream &is, SnapshotTimeQueue &tq);
    };

    
  }
}

#endif
