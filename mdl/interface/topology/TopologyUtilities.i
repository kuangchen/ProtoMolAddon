%module TopologyUtilities
%{
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include "Topo.h"
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/type/Vector3DBlock.h>
%include <protomol/topology/TopologyUtilities.h>
%include "Topo.h"
