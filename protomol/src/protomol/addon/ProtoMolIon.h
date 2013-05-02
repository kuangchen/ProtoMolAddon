#ifndef _PROTOMOL_ION_H
#define _PROTOMOL_ION_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/NeutralAtom.h>

using namespace ProtoMol;
using namespace std;

namespace ProtoMolAddon {

  class ProtoMolIon {
  private:
    double mass;
    double charge;
    Vector3D position;
    Vector3D velocity;
    unsigned int atomId;

  public:
    ProtoMolIon(const ProtoMolApp *app, unsigned int i);
    void UpdateProtoMol(ProtoMolApp *app);

    friend class NeutralAtom;
  };
}

#endif
