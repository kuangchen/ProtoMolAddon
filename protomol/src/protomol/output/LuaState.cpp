#include <protomol/output/LuaState.h>
#include <string>
#include <vector>

using std::vector;
using std::string;
using LuaState::LuaState;

template <> float LuaState::LuaState::lua_get<float>() {
  return lua_tonumber(L, -1);
}

template <> double LuaState::LuaState::lua_get<double>() {
  return lua_tonumber(L, -1);
}

template <> bool LuaState::LuaState::lua_get<bool>() {
  return lua_toboolean(L, -1);
}

template <> int LuaState::LuaState::lua_get<int>() {
  return lua_tointeger(L, -1);
}

template <> string LuaState::LuaState::lua_get<string>() {
  return string(lua_tostring(L, -1));
}

template <> std::vector<double> LuaState::LuaState::lua_get<std::vector<double> > () {
  vector<double> value;

  if (lua_type(L, -1) != LUA_TTABLE)
      luaL_error(L, "Not an array!");

  int len = lua_objlen(L, -1);
  value.resize(len);

  for (int i = 1; i <= len; i++){
    lua_rawgeti(L, -1, i);
    value[i-1] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  return value;
}
