#ifndef _LUA_CONFIG_READER_H
#define _LUA_CONFIG_READER_H

#ifdef ARCHLINUX

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#else
extern "C" {
#include <lua5.2/lua.h>
#include <lua5.2/lauxlib.h>
#include <lua5.2/lualib.h>
}

#endif
#include <sstream>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <exception>

using namespace std;

namespace ProtoMolAddon {

  namespace Lua {

    class LuaConfigReaderException : public std::exception {
    private:
      string msg;

    public:
      LuaConfigReaderException(const string& message) : msg(message) {};
      ~LuaConfigReaderException() throw() {};

      virtual const char* what() const throw() {
	return msg.c_str();
      }

    };


    class LuaConfigReader {
    private:
      lua_State *L;
      template<typename T> T GetValueLowLevel(const char *varname);
      
    public:
      LuaConfigReader(const string& filename);
      LuaConfigReader();
      ~LuaConfigReader();
      
      template<typename T> T GetValue(const char *varname) {
	char temp[1000];
	memset(temp, 0, sizeof(temp));
	int i=0;
	int j=0;
	int n=0;
	while (varname[i] != '\0') {
	  char c = varname[i];
	  if (c == '.') {

	    n==0 ? lua_getglobal(L, temp) : lua_getfield(L, -1, temp);
	    ++n;

	    memset(temp, 0, sizeof(temp));
	    j = 0;
 
	    if (lua_isnil(L, -1)) {
	      ostringstream s;
	      s << "Cannot find field " << varname;
	      throw LuaConfigReaderException(s.str());

	      return T();
	    }
	  }

	  else 
	    temp[j++] = c;

	  ++i;
	}

	n==0 ? lua_getglobal(L, temp) : lua_getfield(L, -1, temp);
	T r = GetValueLowLevel<T>(varname);
	lua_pop(L, n+1);
	return r;
      }

    };

    template <> double LuaConfigReader::GetValueLowLevel<double>(const char *varname);
    template <> bool LuaConfigReader::GetValueLowLevel<bool>(const char *varname);
    template <> int LuaConfigReader::GetValueLowLevel<int>(const char *varname);
    template <> float LuaConfigReader::GetValueLowLevel<float>(const char *varname);
    template <> string LuaConfigReader::GetValueLowLevel<string>(const char *varname);
    template <> vector<double> LuaConfigReader::GetValueLowLevel<vector<double> > (const char *varname);
    template <> vector<string> LuaConfigReader::GetValueLowLevel<vector<string> > (const char *varname);
  }
}

#endif
