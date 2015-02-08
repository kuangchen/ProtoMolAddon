#ifndef _SI_ATOM_PROXY_ARRAY_H
#define _SI_ATOM_PROXY_ARRAY_H

#include <protomol/addon/util/SIAtomProxy.h>
#include <vector>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Util {

    using namespace ProtoMol;
    
    class SIAtomProxyArray : public std::vector<SIAtomProxy> {
    public:
      SIAtomProxyArray(ProtoMol::ProtoMolApp* app);
    };


  }
}

#endif
