#include <protomol/addon/LuaConfigReader.h>
#include <string>
#include <vector>

using namespace ProtoMolAddon::Lua;

LuaConfigReader::LuaConfigReader(const string &filename):
  L(luaL_newstate()) {
  
  luaL_openlibs(L);
  if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
    throw LuaConfigReaderException("Cannot open lua file " + filename);
}

LuaConfigReader::LuaConfigReader(): L(NULL) {}

LuaConfigReader::~LuaConfigReader() {
  //if(L) lua_close(L);
  cout << "Lua Closed\n";
}


template <> float LuaConfigReader::GetValueLowLevel<float>(const char *varname) {
  
  if (!lua_isnumber(L, -1)) {
    ostringstream s;
    s << varname << " is not a float"; 
    throw LuaConfigReaderException(s.str());
  }

  return static_cast<float>(lua_tonumber(L, -1));
}


template <> double LuaConfigReader::GetValueLowLevel<double>(const char *varname) {
  if (!lua_isnumber(L, -1)) {
    ostringstream s;
    s << varname << " is not a double"; 
    throw LuaConfigReaderException(s.str());
  }

  return static_cast<double>(lua_tonumber(L, -1));
}

template <> int LuaConfigReader::GetValueLowLevel<int>(const char* varname) {
  if (!lua_isnumber(L, -1)) {
    ostringstream s;
    s << varname << " is not an int"; 
    throw LuaConfigReaderException(s.str());
  }

  return static_cast<int>(lua_tointeger(L, -1));
}


template <> bool LuaConfigReader::GetValueLowLevel<bool>(const char *varname) {
  if (!lua_isboolean(L, -1)) {
    ostringstream s;
    s << varname << " is not a boolean"; 
    throw LuaConfigReaderException(s.str());
  }
  return static_cast<bool>(lua_toboolean(L, -1));
}

template <> string LuaConfigReader::GetValueLowLevel<string>(const char *varname) {
  if (!lua_isstring(L, -1)) {
    ostringstream s;
    s << varname << " is not a string"; 
    throw LuaConfigReaderException(s.str());
  }
  return string(lua_tostring(L, -1));
}

template <> vector<double> LuaConfigReader::GetValueLowLevel<vector<double> > (const char *varname) {
  if (!lua_istable(L, -1)) {
    ostringstream s;
    s << varname << " is not a table"; 
    throw LuaConfigReaderException(s.str());
  }

  vector<double> value(0);

  //int len = lua_rawlen(L, -1);
  int len = lua_rawlen(L, -1);
  for (int i = 1; i <= len; i++){
    lua_rawgeti(L, -1, i);

    if (lua_isnumber(L, -1))
      value.push_back(lua_tonumber(L, -1));
    
    lua_pop(L, 1);
  }

  return value;
}

template <> vector<string> LuaConfigReader::GetValueLowLevel<vector<string> > (const char *varname) {
  if (!lua_istable(L, -1)) {
    ostringstream s;
    s << varname << " is not a table"; 
    throw LuaConfigReaderException(s.str());
  }

  vector<string> value(0);

  int len = lua_rawlen(L, -1);
  //int len = lua_objlen(L, -1);
  for (int i = 1; i <= len; i++){
    lua_rawgeti(L, -1, i);

    if (!lua_isstring(L, -1)) {
      ostringstream s;
      s << varname << " is not a string"; 
      throw LuaConfigReaderException(s.str());
    }
     
    value.push_back(string(lua_tostring(L, -1)));
    lua_pop(L, 1);
  }

  return value;
}

