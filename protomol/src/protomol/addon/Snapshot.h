#ifndef __SNAPSHOT_H
#define __SNAPSHOT_H

#include <H5Cpp.h>
#include <H5File.h>
#include <iostream>
#include <string>
#include <protomol/ProtoMolApp.h>

using namespace H5;
using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {

  namespace IonSnapshot {

    typedef struct SnapshotConfig {
      double start;
      int fps;
      int numFrame;
      int numAtom;
      int numFrameperMM;
    } SnapshotConfig;

    class Snapshot {
    private:
      string filename;
      int nextFrame;
      H5File file;

      SnapshotConfig config;

      double start;
      int numFrame;
      int numAtom;

    public:
      Snapshot(const string& filename, const SnapshotConfig& config);

      ~Snapshot();
    
      bool IsActive(double time);
      void AddHeader();
      void AddFrame(double time, double *data);
    };

  }
}

#endif
