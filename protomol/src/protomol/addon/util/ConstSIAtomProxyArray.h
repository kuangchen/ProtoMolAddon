#ifndef _CONST_SI_ATOM_PROXY_ARRAY_H
#define _CONST_SI_ATOM_PROXY_ARRAY_H

#include <protomol/addon/util/ConstSIAtomProxy.h>
#include <vector>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Util {

    using namespace ProtoMol;
    
    class ConstSIAtomProxyArray : public std::vector<ConstSIAtomProxy> {
    public:
      ConstSIAtomProxyArray(const ProtoMolApp* app);
    };


  }
}

#endif
