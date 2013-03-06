%module ImproperForce
%{
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/force/bonded/ImproperSystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/bonded/MTorsionSystemForce.h>
%include <protomol/force/bonded/ImproperSystemForce.h>
%template(MTSF_Periodic) ProtoMol::MTorsionSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(MTSF_Vacuum) ProtoMol::MTorsionSystemForce <ProtoMol::VacuumBoundaryConditions >;
%template(ISF_Periodic) ProtoMol::ImproperSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(ISF_Vacuum) ProtoMol::ImproperSystemForce <ProtoMol::VacuumBoundaryConditions >;