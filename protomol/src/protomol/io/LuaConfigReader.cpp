#include "LuaConfigReader.h"

using namespace std;

namespace Util
{
  void ExtractInteger(lua_State *L, int* value, const string& name, const string& error_msg)
  {
    lua_getglobal(L, name.c_str());
    if (lua_type(L, -1) != LUA_TNUMBER)
      luaL_error(L, error_msg.c_str());

    *(value) = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }

  void ExtractDouble(lua_State *L, double* value, const string& name, const string& error_msg)
  {
    lua_getglobal(L, name.c_str());
    if (lua_type(L, -1) != LUA_TNUMBER)
      luaL_error(L, error_msg.c_str());

    *(value) = (double)lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  void ExtractString(lua_State *L, string& value, const string& name, const string& error_msg)
  {
    lua_getglobal(L, name.c_str());
    if (lua_type(L, -1) != LUA_TSTRING)
      luaL_error(L, error_msg.c_str());

    value.assign(lua_tostring(L, -1));
    lua_pop(L, 1);
  }

  void ExtractArray(lua_State *L, vector<double>& value, const string& name, const string& error_msg)
  {
    lua_getglobal(L, name.c_str());
    if (lua_type(L, -1) != LUA_TTABLE)
      luaL_error(L, error_msg.c_str());

    int len = lua_objlen(L, -1);
    value.resize(len);

    for (int i = 1; i <= len; i++){
      lua_rawgeti(L, -1, i);
      value[i-1] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }

    lua_pop(L, 1);
  }
}
