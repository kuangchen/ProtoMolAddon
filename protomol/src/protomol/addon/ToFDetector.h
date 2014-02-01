#ifndef __TOF_DETECTOR_H
#define __TOF_DETECTOR_H

#include <protomol/ProtoMolApp.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <string>
#include <vector>
#include <memory>

using std::shared_ptr;
using std::string;
using std::vector;
using std::ostream;
using namespace ProtoMol;

namespace ProtoMolAddon {

  class ToFDetector {
  public:
    class AtomRecord {
    public:
      const string &name;
      const Vector3D &pos;
      const Vector3D &vel;
      double time;
      bool hit;

    private:
      AtomRecord(const AtomRecord &r);
      const AtomRecord& operator= (const AtomRecord &r);

    public:
      AtomRecord(const string &name, const Vector3D &pos, const Vector3D &vel, double time=0);

      bool Hit() const { return hit; }
      void SetHitTime(double t) { hit = true; time = t; }

      friend ostream& operator<< (ostream &os, const AtomRecord &r);
      friend class ToFDetector;
    };


  private:
    vector< shared_ptr<AtomRecord> > record_list;
    double detector_pos;

  public:
    ToFDetector(double detector_pos=0);

    void Initialize(const ProtoMolApp *app);
    void UpdateRecord(const ProtoMolApp *app);
    friend ostream& operator<< (ostream &os, const ToFDetector &detector);
    friend class OutputToF;
  };

}

#endif
 
