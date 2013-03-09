extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include <protomol/addon/LinearPaulTrap.h>
#include <iostream>

using namespace ProtoMolAddon::LinearPaulTrap;
using namespace std;
  
Lqt::Lqt() {}

Lqt::Lqt(LuaState& L) {

  m = L.get<double>("trap.m");
  r0 = L.get<double>("trap.r0");
  z0 = L.get<double>("trap.z0");
  v_rf = L.get<double>("trap.v_rf");
  v_ec = L.get<double>("trap.v_ec");
  eta = L.get<double>("trap.eta");
  omega = L.get<double>("trap.omega");

  double freq_sq = omega * omega * 4 * M_PI * M_PI;

  double qq = -4 * v_rf / (r0 * r0 * freq_sq) / (m * 1.0364e-8);
  double aa = -4 * v_ec * eta / (z0 * z0 * freq_sq) / (m * 1.0364e-8) * 4;

  q.push_back(qq);
  q.push_back(-qq);
  q.push_back(0);

  a.push_back(aa);
  a.push_back(aa);
  a.push_back(-2*aa);
    
  for (int i=0; i<3; i++)
    mf.push_back(MathieuFunction(q[i], a[i], 15));
		   
}

void Lqt::GetEnergy(const Vector3D& pos, const Vector3D& vel, double time, vector<double>& totEnergy, vector<double>& secEnergy) {
  double scale = M_PI * omega;
  MathieuMiscVar mmv;
  double w;
  vector<double> r;
  r.resize(4);
  totEnergy.resize(0);
  secEnergy.resize(0);
  double e;
	
  for(int i=0; i<3; i++) {
    mf[i].GetMiscVar(&mmv);
    w = mmv.wronskian;
    mf[i].Evaluate(time*scale, r);
    e = ((r[2] * r[2] + r[3] * r[3]) * pos[i] * pos[i] +  
	 (r[0] * r[0] + r[1] * r[1]) * vel[i] * vel[i] / scale / scale  +
	 (r[0] * r[2] + r[1] * r[3]) * pos[i] * vel[i] * (-2) / scale) / w / w * mmv.cp_sp_sq_sum_ave / 2.0 * m * 1.66e-27 / 1.38e-23 * scale * scale;
    totEnergy.push_back(e);
    secEnergy.push_back(e*mmv.eta);
  }
}

void Lqt::GetForce(const Vector3D& pos, double time, Vector3D& force)
{
  double a_rf = 2 * 1.60217646e-19 / r0 / r0 * v_rf * cos(2*M_PI*omega*time);
  double a_dc = eta * 1.60217646e-19 / z0 / z0 * 4 * v_ec;

  force[0] =   -a_rf * pos[0] + a_dc * pos[0];
  force[1] =    a_rf * pos[1] + a_dc * pos[1];
  force[2] = -2*a_dc * pos[2];
}

vector<double> Lqt::GetSecularFrequency() {
  vector<double> w(3, 0);

  for (int i=0; i<3; i++)
    w[i] = omega / 2 * mf[i].GetMu();

  return w;
}

double Lqt::GetFrequency() {
  return omega;
}
