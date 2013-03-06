%module AngleForce
%{
#include <protomol/topology/Angle.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/force/bonded/AngleSystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/bonded/AngleSystemForce.h>
%template(ASF_Periodic) ProtoMol::AngleSystemForce<ProtoMol::PeriodicBoundaryConditions>;
%template(ASF_Vacuum) ProtoMol::AngleSystemForce<ProtoMol::VacuumBoundaryConditions>;