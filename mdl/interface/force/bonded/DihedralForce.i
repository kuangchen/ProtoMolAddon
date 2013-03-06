%module DihedralForce
%{
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/force/bonded/DihedralSystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}


%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/bonded/MTorsionSystemForce.h>
%include <protomol/force/bonded/DihedralSystemForce.h>
%template(MTSF_Periodic) ProtoMol::MTorsionSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(MTSF_Vacuum) ProtoMol::MTorsionSystemForce <ProtoMol::VacuumBoundaryConditions >;
%template(DSF_Periodic) ProtoMol::DihedralSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(DSF_Vacuum) ProtoMol::DihedralSystemForce <ProtoMol::VacuumBoundaryConditions >;