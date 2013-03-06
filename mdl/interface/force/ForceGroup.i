%module ForceGroup
%{
#include <protomol/force/ForceGroup.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/force/system/SystemForce.h>
using namespace ProtoMol;
%}

%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/ForceGroup.h>

%extend ProtoMol::ForceGroup {
    void clear() {
      self->getForces().clear();
    }
}