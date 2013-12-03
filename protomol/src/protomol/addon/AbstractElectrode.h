#ifndef _ABSTRACT_ELECTRODE_H
#define _ABSTRACT_ELECTRODE_H

#include "boost/multi_array.hpp"
#include <protomol/type/Vector3D.h>
#include <vector>
#include <string>
#include <bitset>
#include <array>

using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {

  class AbstractElectrode {
  public:
    string label;

  public:
    AbstractElectrode(const string& label): label(label) {};
    virtual ~AbstractElectrode() {};
virtual double GetNNPotential(const Vector3D& pos, const boost::array<int, 3>& offset=boost::array<int, 3>()) const = 0;
    virtual double GetNNVoltage(double time, int offset=0) const = 0;
    virtual const Vector3D& GetDx() = 0;
    virtual void GetFraction(const Vector3D &pos, boost::array<double, 3>& f) = 0;
    virtual void DumpInfo(ostream& os) = 0;
    
  };
}

#endif
