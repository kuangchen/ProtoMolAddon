#ifndef _TOF_H
#define _TOF_H

#include <protomol/type/Vector3D.h>
#include <protomol/base/PMConstants.h>
#include <vector>
#include <string>
#include <array>
#include <protomol/addon/AbstractElectrode.h>
#include <protomol/addon/Electrode.h>

using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {
  
  class ToF {
  private:
    vector< Electrode > elct;

  public:
    ToF();
    ToF(const string &def);
    ~ToF();
    double GetTotalRealTimePotential(const Vector3D& pos, double t, const boost::array<int, 3>& offset);
    double GetTotalRealTimeInterpolatedPotential(const Vector3D& pos, double t);
    void GetForce(double charge, const Vector3D &pos, double t, Vector3D& force);
    void DumpElectrodeInfo(ostream& os);

    friend ostream& operator<< (ostream& os, const ToF& tof);

    static void Test();
  };

}

#endif
