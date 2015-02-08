#include <protomol/addon/util/ConstSIAtomProxyArray.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMolAddon {
  namespace Util {

    ConstSIAtomProxyArray::ConstSIAtomProxyArray(const ProtoMolApp *app) {
      for (unsigned int i=0; i<app->positions.size(); i++)
	push_back(ConstSIAtomProxy(app, i));
    }

  }
}
