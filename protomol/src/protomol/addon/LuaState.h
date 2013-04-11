#ifndef _LUA_STATE_H
#define _LUA_STATE_H

#ifdef ARCHLINUX

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#else
extern "C" {
#include <lua5.1/lua.h>
#include <lua5.1/lauxlib.h>
#include <lua5.1/lualib.h>
}


#endif

#include <string>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

namespace ProtoMolAddon {

  namespace Lua {

    class LuaStateError {
    private:
      string msg;
    public:
      LuaStateError(const string& message);
      string GetMessage() const;
    };
    
    class LuaState {
    private:
      lua_State *L;
 
    public:
      LuaState(const string& filename) {
	L = luaL_newstate();
	luaL_openlibs(L);
	if (luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
	  luaL_error(L, "cannot run configuration file: %s",
		     lua_tostring(L, -1));
      }
    
      LuaState() :
	L(NULL) {}

      ~LuaState() {
	if (L) 
	  lua_close(L);
      }

      template<typename T> T get(const char *varname) {
	char temp[64];
	memset(temp, 0, sizeof(temp));
	int i=0;
	int j=0;
	int n=0;
	while (varname[i] != '\0') {
	  char c = varname[i];
	  if (c == '.') {

	    if (n == 0)  // This snippet push field onto the stack
	      lua_getglobal(L, temp);
	    
	    else
	      lua_getfield(L, -1, temp);

	    ++n;

	    memset(temp, 0, sizeof(temp));
	    j = 0;
 
	    if (lua_isnil(L, -1)) {
	      throw LuaStateError( "Error reading field " + string(varname) );
	      return T();
	    }
	  }

	  else {
	    temp[j] = c;
	    ++j;
	  }
	  ++i;
	}

	if (n == 0)
	  lua_getglobal(L, temp);
	else
	  lua_getfield(L, -1, temp);
      
	return lua_get<T>();
      }

      // Generic get
      template<typename T> T lua_get() {
	return 0;
      }
    };

    template <> double LuaState::lua_get<double>();
    template <> bool LuaState::lua_get<bool>();
    template <> int LuaState::lua_get<int>();
    template <> float LuaState::lua_get<float>();
    template <> string LuaState::lua_get<string>();
    template <> vector<double> LuaState::lua_get<vector<double> > ();
  }

  

}
#endif
