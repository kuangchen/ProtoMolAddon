#ifndef _TOF_H
#define _TOF_H

#include <protomol/type/Vector3D.h>
#include <protomol/base/PMConstants.h>
#include <vector>
#include <string>
#include <array>
#include <protomol/addon/AbstractElectrode.h>

using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {
  

  class ToF {
  private:
    vector< vector<AbstractElectrode*> > electrode;

  public:
    ToF();
    ToF(const string &def);
    ~ToF();
    double GetTotalPotential(const Vector3D& pos, double t,  const array<int, 3>& offset=array<int, 3>());
    void GetForce(double charge, const Vector3D &pos, double t, Vector3D& force);
    void DumpElectrodeInfo(ostream& os);

    static void Test();
  };

}

#endif
