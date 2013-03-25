#include <protomol/addon/LuaState.h>
#include <string>
#include <vector>

using namespace ProtoMolAddon::Lua;

template <> float LuaState::lua_get<float>() {
  return lua_tonumber(L, -1);
}

template <> double LuaState::lua_get<double>() {
  return lua_tonumber(L, -1);
}

template <> bool LuaState::lua_get<bool>() {
  return lua_toboolean(L, -1);
}

template <> int LuaState::lua_get<int>() {
  return lua_tointeger(L, -1);
}

template <> string LuaState::lua_get<string>() {
  return string(lua_tostring(L, -1));
}

template <> vector<double> LuaState::lua_get<vector<double> > () {
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

LuaStateError::LuaStateError(const string& message) : msg(message) {
  std::cout << message << "\n";
}

string LuaStateError::GetMessage() const {
  return msg;
}
