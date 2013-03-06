#ifndef __SNAPSHOT_H
#define __SNAPSHOT_H

#include <H5Cpp.h>
#include <H5File.h>
#include <iostream>
#include <string>
#include <protomol/ProtoMolApp.h>


using namespace H5;
using namespace std;

namespace ProtoMol {

  typedef struct SnapshotConfig {
    Real start;
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

    Real start;
    int numFrame;
    int numAtom;

  public:
    Snapshot(const string& filename, const SnapshotConfig& config);

    ~Snapshot();
    
    bool IsActive(Real time);
    void AddHeader();
    void AddFrame(Real time, Real *data);
  };
}

#endif
