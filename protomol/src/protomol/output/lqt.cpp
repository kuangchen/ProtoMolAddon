extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include "lqt.h"
#include <iostream>

using namespace std;

namespace LQT
{

  typedef struct
  {
    double *p_mass;
    LQT* p_lqt;
  }param_struct;
    


  double LQT::ExtractParameter(lua_State *L, const std::string& name, const std::string& error_msg)
  {
    lua_getglobal(L, name.c_str());
     
    if (!lua_isnumber(L, -1)) 
      luaL_error(L, error_msg.c_str());

    double r = (double)lua_tonumber(L, -1);
    lua_pop(L, 1);
    return r;
  }

  LQT::LQT()
  {}

  LQT::LQT(double radius, double freq, double volt_rf, double volt_ec, double phase, double length, double kappa):
    radius(radius), freq(freq), volt_rf(volt_rf), volt_ec(volt_ec), phase(phase), length(length), kappa(kappa)
  {
    
  }

  LQT::LQT(const std::string& def_filename)
  {
    std:: cout << "opening " << def_filename <<"\n";
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);
    
    if (luaL_loadfile(L, def_filename.c_str()) || lua_pcall(L, 0, 0, 0))
      luaL_error(L, "cannot run configuration file: %s",
		   lua_tostring(L, -1));

    radius  = ExtractParameter(L, "radius",  "radius should be a number\n" );
    freq    = ExtractParameter(L, "freq",    "freq should be a number\n"   );
    volt_rf = ExtractParameter(L, "volt_rf", "volt_rf should be a number\n");
    volt_ec = ExtractParameter(L, "volt_ec", "volt_ec should be a number\n");
    phase   = ExtractParameter(L, "phase",   "phase should be a number\n"  );
    length  = ExtractParameter(L, "length",  "length should be a number\n" );
    kappa   = ExtractParameter(L, "kappa",   "kappa should be a number\n"  ); 
    
    lua_close(L);

    
  }

  void LQT::GetForce(double mass_charge_ratio, const Vector3D& position, const Vector3D& velocity, double time, Vector3D& force)
  {
    double prefactor_rf = 2 * 1.60217646e-19 / radius / radius * volt_rf * cos(2*3.14159265*freq*time + phase);
    double prefactor_dc = kappa * 1.60217646e-19 / length / length * 4 * volt_ec;

    force[0] = -prefactor_rf * position[0] + prefactor_dc * position[0];
    force[1] =  prefactor_rf * position[1] + prefactor_dc * position[1];
    force[2] = -2*prefactor_dc * position[2];
  }
}
