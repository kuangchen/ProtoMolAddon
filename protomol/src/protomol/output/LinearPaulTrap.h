#ifndef _LQT2_H
#define _LQT2_H

#include <vector>
#include <protomol/output/Mathieu.h>
#include <protomol/type/Vector3D.h>
#include <protomol/output/LuaState.h>

using namespace LuaState;
using namespace ProtoMol;
using namespace Mathieu;
using namespace std;

namespace LinearPaulTrap {

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
    Lqt(LuaState::LuaState& L);
    Lqt();
    ~Lqt() {};

    int Initialize(LuaState::LuaState& L);
    double GetFrequency();
    void GetForce(const Vector3D& pos, double time, Vector3D& force);
    void GetEnergy(const Vector3D& pos, const Vector3D& vel, double time, vector<double>& totEnergy, vector<double>& secEnergy);
    vector<double> GetSecularFrequency();
  };
}

#endif
