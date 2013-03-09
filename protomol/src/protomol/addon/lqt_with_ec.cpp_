extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include "lqt_with_ec.h"
#include <iostream>
#include <protomol/io/LuaConfigReader.h>

using namespace std;
using namespace Util;

namespace LQTWithEC
{
  LQTWithEC::LQTWithEC()
  {}

  LQTWithEC::LQTWithEC(const std::string& def_filename)
  {
    std:: cout << "opening " << def_filename <<"\n";
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);
    
    if (luaL_loadfile(L, def_filename.c_str()) || lua_pcall(L, 0, 0, 0))
      luaL_error(L, "cannot run configuration file: %s",
		   lua_tostring(L, -1));

    ExtractDouble(L, &radius, "radius",   "radius should be a number\n" );
    ExtractDouble(L, &freq,   "freq",     "freq should be a number\n"   );
    ExtractDouble(L, &volt_rf,"volt_rf",  "volt_rf should be a number\n");
    ExtractDouble(L, &volt_ec,"volt_ec",  "volt_ec should be a number\n");
    ExtractDouble(L, &phase,  "phase",    "phase should be a number\n"  );
    ExtractDouble(L, &length, "length",   "length should be a number\n" );
    ExtractArray(L, kappa_rf,"kappa_rf", "error reading kappa_rf\n"    );
    ExtractArray(L, kappa_ec,"kappa_ec", "error reading kappa_ec\n"    );

    lua_close(L);

    r0_square = radius * radius;
    z0_square = length * length / 4;
  }

  void LQTWithEC::GetForce(const Vector3D& position, const Vector3D& velocity, double time, Vector3D& force){

    double volt_rf_t = volt_rf * cos(2*3.14159265*freq*time + phase);
    coeff[0] = volt_rf_t * kappa_rf[0] + volt_ec * kappa_ec[0];
    coeff[1] = volt_rf_t * kappa_rf[1] + volt_ec * kappa_ec[1];

    vector<double> gradient[2];
    gradient[0].resize(3);
    gradient[0][0] = -2 * position[0];
    gradient[0][1] =  2 * position[1];
    gradient[0][2] =  0;

    gradient[1].resize(3);
    gradient[1][0] =      position[0];
    gradient[1][1] =      position[1];
    gradient[1][2] = -2 * position[2];
    
    for (unsigned int i = 0 ; i < 3 ; i++)
      force[i] = (coeff[0] * gradient[0][i] / r0_square + coeff[1] * gradient[1][i] / z0_square) * 1.60217649e-19;
  }
}
