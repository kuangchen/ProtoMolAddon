#ifndef __DUMMY_ELECTRODE_H
#define __DUMMY_ELECTRODE_H

#include <protomol/type/Vector3D.h>

// extern "C" {
// #include <lua5.2/lua.h>
// #include <lua5.2/lualib.h>
// #include <lua5.2/lauxlib.h>
// }

// #include <luabind/luabind.hpp>
// #include <luabind/function.hpp>

namespace ProtoMolAddon {
  namespace SegmentedTrap {
    
    using namespace ProtoMol;

    class DummyElectrode {
    private:

    public:
      DummyElectrode();
      double GetPotential(const ProtoMol::Vector3D &pos, double t) const;

    };
  }
}


#endif
