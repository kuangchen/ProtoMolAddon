#ifndef _LQT2_H
#define _LQT2_H

#include <vector>
#include <protomol/addon/Mathieu.h>
#include <protomol/addon/LuaState.h>
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
      Lqt(LuaState& L);
      Lqt();
      ~Lqt() {};

      double GetFrequency();
      void GetForce(const Vector3D& pos, double time, Vector3D& force);
      void GetEnergy(const Vector3D& pos, const Vector3D& vel, double time, vector<double>& totEnergy, vector<double>& secEnergy);
      vector<double> GetSecularFrequency();
    };
  }
}
#endif
