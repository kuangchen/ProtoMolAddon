#ifndef _LQT_With_EC_H
#define _LQT_With_EC_H

#include <string>
#include <protomol/type/Vector3D.h>

namespace LQTWithEC
{
  using namespace ProtoMol;

  class LQTWithEC
  {
  private:
    double radius;
    double freq;
    double volt_rf;
    double volt_ec;
    double phase;
    double length;
    vector<double> kappa_rf, kappa_ec;

    double coeff[2];
    double r0_square;
    double z0_square;

  public:
    LQTWithEC(const std::string& def_filename);
    LQTWithEC();

  public:
    void GetForce(const Vector3D& position, const Vector3D& velocity, double time, Vector3D& force);

  };
}

#endif
