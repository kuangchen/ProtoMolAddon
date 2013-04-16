#ifndef _LQT2_H
#define _LQT2_H

#include <vector>
#include <protomol/addon/Mathieu.h>
#include <protomol/addon/LuaConfigReader.h>
#include <protomol/type/Vector3D.h>

using namespace ProtoMolAddon::Lua;
using namespace ProtoMolAddon::Mathieu;
using namespace ProtoMol;
using namespace std;

namespace ProtoMolAddon {

  namespace LinearPaulTrap {
    using namespace std;

    class Lqt {
    private:
      double m;
      double r0;
      double z0;
      double v_rf;
      double v_ec;
      double eta;
      double omega;

      vector<double> q;
      vector<double> a;
      vector<MathieuFunction> mf;

    public:
      Lqt(LuaConfigReader& reader);
      Lqt();
      ~Lqt() {};

      double GetFrequency();
      Vector3D GetForce(double charge, const Vector3D& pos, double time);
      void GetEnergy(const Vector3D& pos, const Vector3D& vel, double time, vector<double>& totEnergy, vector<double>& secEnergy);
      vector<double> GetSecularFrequency();
    };
  }
}
#endif
