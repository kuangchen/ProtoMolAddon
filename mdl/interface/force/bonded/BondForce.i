%module BondForce
%{
#include <protomol/topology/Bond.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/force/bonded/BondSystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}
%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/bonded/BondSystemForce.h>


%template(BSF_Periodic) ProtoMol::BondSystemForce<ProtoMol::PeriodicBoundaryConditions >;
%template(BSF_Vacuum) ProtoMol::BondSystemForce<ProtoMol::VacuumBoundaryConditions >;
