#ifndef __TOF_DETECTOR_H
#define __TOF_DETECTOR_H

#include <protomol/ProtoMolApp.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMol;

namespace ProtoMolAddon {

  class ToFDetector {
  private:
    vector<bool> hit;
    vector<double> hit_time;
    vector<Vector3D> hit_position;
    vector<Vector3D> hit_velocity;
    double detector_pos;

  public:
    ToFDetector(double detector_pos=0);
    ~ToFDetector();
    void Initialize(const ProtoMolApp* app);
    void Update(const ProtoMolApp* app);
    void Finalize();
    friend ostream& operator<<(ostream& os, ToFDetector& tof);
    friend class OutputToF;
  };

}

#endif
 
