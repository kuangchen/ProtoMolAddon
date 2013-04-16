#include <protomol/addon/NeutralAtom.h>
#include <protomol/addon/NeutralAtomLookupTable.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/Util.h>
#include <cmath>

using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Util;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Constant;

NeutralAtom::NeutralAtom(double m, double pl, double t, double d, const Vector3D &pos, const Vector3D &vel) :
  mass(m), polarizability(pl), temperature(t), density(d), position(pos), velocity(vel),
  rd(),
  uniform_dist(0, 1),
  mag_dist(),
  polar_angle_dist(),
  azimuthal_angle_dist(0, 2*M_PI) {

   if (mass>1e-50) 
     sigma = sqrt(SI::BOLTZMANN*temperature/mass);
 
   LoadLookupTable();
   LoadCache();
}
  

NeutralAtom::NeutralAtom(LuaConfigReader &reader):
  mass(reader.GetValue<double>("neutral.mass") * SI::AMU),
  polarizability(reader.GetValue<double>("neutral.polarizability") * 1e-30),
  temperature(reader.GetValue<double>("neutral.temperature")),
  density(reader.GetValue<double>("neutral.density")),
  position(Vector3D()), 
  velocity(Vector3D()) 
{
  if (mass>1e-50) 
    sigma = sqrt(SI::BOLTZMANN*temperature/mass);
 
   LoadLookupTable();
   LoadCache();

}
  
void NeutralAtom::CollideAll(ProtoMolApp *app, double dt) {
  for (unsigned int i=0; i<app->velocities.size(); i++) {
    double q = app->topology->atoms[i].scaledCharge / SQRTCOULOMBCONSTANT * SI::ELECTRON_CHARGE;
    double m = app->topology->atoms[i].scaledMass * SI::AMU;
    Vector3D v = app->velocities[i] * VELOCITY_CONV;

    double rate = GetAverageCollisionRate(m, q, v);
    double dice = uniform_dist(rd);

    if (rate * dt > dice) {
      //      std::cout << "Collide " << i << "\n";

      velocity = SampleVelocity(v);
      velocity = Rotate(v, velocity);
      Collide(m, v);
      app->velocities[i] = v / VELOCITY_CONV;
    }
  }
}

double NeutralAtom::GetCollisionRate(double m, double q, const Vector3D &v) {
  double C4 = polarizability * q * q / (4*M_PI*EPSILON_0);
  double prefactor = M_PI * pow(pow(C4/HBAR, 2)*sigma*2, 1.0/3) * (1+M_PI*M_PI/16) * density;

  double v_norm = (v-velocity).norm() / sigma;  
  return prefactor * pow(v_norm, 1.0/3);
}

double NeutralAtom::GetLangevinRate(double m, double q, const Vector3D &v) {
  double C4 = polarizability * q * q / (4*M_PI*EPSILON_0);
  double mu = m * mass / (m+mass);
  return 2 * M_PI * sqrt(C4 / mu) * density;
}

double NeutralAtom::GetAverageCollisionRate(double m, double q, const Vector3D &v) {
  double C4 = polarizability * q * q / (4*M_PI*EPSILON_0);
  double prefactor = M_PI * pow(pow(C4/HBAR, 2)*sigma*2, 1.0/3) * (1+M_PI*M_PI/16) * density;

  double v_norm = v.norm() / sigma;  
  if (v_norm > 5) 
    return prefactor * pow(v_norm, 1.0/3);

  else {
    std::map<double,double>::iterator it1 = lookup_table1.upper_bound(v_norm);
    std::map<double,double>::iterator it2 = it1;
    --it1;    

    // std::cout << "\nlower bound = " << it1->first;
    // std::cout << "\nupper bound = " << it2->first;

    double f = (v_norm - it1->first) / (it2->first - it1->first);

    return (f * it2->second + (1-f) * it1->second) * prefactor;
  }
}

void NeutralAtom::LoadLookupTable() {
  int i=0;
  double vi, integral;

  while (1) {
    vi = table1[i][0];
    integral= table1[i][1];
    
    if (vi<0 || integral<0) break;
    std::cout << "Adding key = " << vi << "\t value = " << integral << "\n";
    lookup_table1[vi] = integral;
    i++;
  }
}

Vector3D NeutralAtom::SampleVelocity(const Vector3D &v) {
  typedef std::piecewise_constant_distribution<double>::param_type param_type;

  double vi = v.norm() / sigma;

  for (double r:r_cache_middle) 
    w_r_cache.push_back(r*exp(-r*r/2)*(pow((r+vi),7.0/3)-pow(fabs(r-vi),7.0/3)));
  
  param_type pr(r_cache.begin(), r_cache.end(), w_r_cache.begin());
  mag_dist.param(pr);
  double r = mag_dist(rd);
  w_r_cache.resize(0);

  for (double theta:theta_cache_middle)
    w_theta_cache.push_back(pow((r*r-2*r*vi*cos(theta)+vi*vi), 1.0/6));
  
  param_type pt(theta_cache.begin(), theta_cache.end(), w_theta_cache.begin());
  polar_angle_dist.param(pt);
  double theta = polar_angle_dist(rd);
  w_theta_cache.resize(0);

  double phi = azimuthal_angle_dist(rd);
  
  return BuildVector(r*sigma, theta, phi);
}


void NeutralAtom::Collide(double m, Vector3D &v) {
  double b1 = mass / (mass + m);
  double b2 = 1-b1;

  Vector3D v_com, v_rel, v_rel_after;
  v_com = v * b1 + velocity * b2;
  v_rel = v - velocity;

  double r = v_rel.norm();

  double theta = M_PI;
  
  double phi = azimuthal_angle_dist(rd);

  v_rel_after = BuildVector(r, theta, phi);
  v_rel_after = Rotate(v_rel, v_rel_after);
  //std::cout << "v_rel_before = " << v_rel << "\t" << "after = " << v_rel_after << "\n";
  v = v_com + v_rel_after * b2;
}


void NeutralAtom::LoadCache() {
  r_cache.resize(101);
  theta_cache.resize(101);

  r_cache_middle.resize(100);
  theta_cache_middle.resize(100);

  for (int i=0; i<100; i++) {
    r_cache[i] = 0.1 * i; 
    r_cache_middle[i] = r_cache[i] + 0.05;

    theta_cache[i] = 0.01 * M_PI * i;
    theta_cache_middle[i] = theta_cache[i] + 0.01 * M_PI * 0.5;
  }

  r_cache[100] = 10;
  theta_cache[100] = M_PI;

  w_r_cache.resize(0);
  w_theta_cache.resize(0);

}

void NeutralAtom::Test() {
  const int ave = 10000;
  Vector3D v(-18, 100, -100);

  double rate = 0;

  for (int i=0; i<ave; i++) {
    if (i%100==0) std::cout << i/100 << "%" << "\n";
    velocity = SampleVelocity(v);
    velocity = Rotate(v, velocity);
    rate += GetCollisionRate(173*SI::AMU, SI::ELECTRON_CHARGE, v);
  }
  
  rate /= ave;
  std::cout << "Rate from Monte Carlo " << rate << "\n";
  std::cout << "Rate from Langevin    " << GetLangevinRate(173*SI::AMU, SI::ELECTRON_CHARGE, v) << "\n";
  std::cout << "Rate from Direct Calc " << GetAverageCollisionRate(173*SI::AMU, SI::ELECTRON_CHARGE, v) << "\n";

}
