#ifndef OUTPUT_ION_SNAPSHOT_H
#define OUTPUT_ION_SNAPSHOT_H

#include <protomol/output/Output.h>
#include <protomol/addon/snapshot/SnapshotManager.h>
#include <protomol/addon/snapshot/HDF5CompressedStorage.h>
#include <protomol/addon/snapshot/EqualSpacingTimeQueue.h>
#include <vector>
#include <string>

using namespace std;
using namespace ProtoMol;
using ProtoMolAddon::Snapshot::SnapshotManager;

namespace ProtoMolAddon {
  namespace Snapshot {

    class OutputIonSnapshot : public Output {
    private:
      string fname;
      SnapshotManager<EqualSpacingTimeQueue, HDF5CompressedStorage> manager;

    public:
      OutputIonSnapshot(const string& filename);
      OutputIonSnapshot();
      ~OutputIonSnapshot();
  
      Output *doMake(const vector<Value> &values) const;
      void doInitialize();
      void doRun(int step);
      void doFinalize(int step);
  
    public:
      static const string keyword;
      void getParameters(vector<Parameter> &parameter) const;
      string getIdNoAlias() const {return keyword;};
    
    };

  }
}

#endif
