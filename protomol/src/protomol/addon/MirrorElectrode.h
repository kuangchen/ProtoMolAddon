#ifndef __MIRROR_ELECTRODE_H
#define __MIRROR_ELECTRODE_H

#include <protomol/addon/AbstractElectrode.h>
#include <protomol/addon/RealElectrode.h>

#include <bitset>
#include <string>

using namespace std;

namespace ProtoMolAddon {

  class MirrorElectrode : public AbstractElectrode {
  private:
    const RealElectrode* parent;
    bitset<3> reflection;

  public:
    MirrorElectrode(const string& label, const RealElectrode* parent, const std::string& reflection);
    ~MirrorElectrode() {};
    double GetNNPotential(const Vector3D& pos, const array<int, 3>& offset) const;
    double GetNNVoltage(double time, int offset) const;
    void GetFraction(const Vector3D &pos, array<double, 3>& f);
    const Vector3D& GetDx();
    void DumpInfo(ostream& os);
  };
}

#endif
