#ifndef _LUA_CONFIG_READER_H
#define _LUA_CONFIG_READER_H

extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}

#include <vector>
#include <string>

namespace Util
{
  using namespace std;
  void ExtractInteger(lua_State *L, int* value, const string& name, const string& error_msg);
  void ExtractDouble(lua_State *L, double* value, const string& name, const string& error_msg);
  void ExtractString(lua_State *L, string& value, const string& name, const string& error_msg);
  void ExtractArray(lua_State *L, vector<double>& value, const string& name, const string& error_msg);
}

#endif
