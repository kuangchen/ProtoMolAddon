#include <protomol/addon/util/SIAtomProxyArray.h>

namespace ProtoMolAddon {
  namespace Util {

    SIAtomProxyArray::SIAtomProxyArray(ProtoMolApp *app) {
      if (app!=NULL)
	for (unsigned int i=0; i<app->positions.size(); i++)
	  push_back(SIAtomProxy(app, i));
    }

  }
}
