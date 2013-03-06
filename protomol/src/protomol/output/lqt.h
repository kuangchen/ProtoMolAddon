#ifndef _LQT_H
#define _LQT_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include <string>
#include <protomol/type/Vector3D.h>

namespace LQT
{

  using namespace ProtoMol;

  class LQT
  {
  public:
    double radius;
    double freq;
    double volt_rf;
    double volt_ec;
    double phase;
    double length;
    double kappa;

    double qx, ax;


  public:
    LQT(double radius, double freq, double volt_rf, double volt_ec, double phase, double length, double kappa);
    LQT(const std::string& def_filename);
    LQT();

    double ExtractParameter(lua_State *L, const std::string& name, const std::string& error_message);
    
  public:
    void GetForce(double mass_charge_ratio, const Vector3D& position, const Vector3D& velocity, double time, Vector3D& force);

  };
}

#endif
