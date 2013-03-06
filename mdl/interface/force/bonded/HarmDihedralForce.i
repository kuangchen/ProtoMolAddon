%module HarmDihedralForce
%{
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/force/bonded/HarmDihedralSystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/config/Value.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}



%include <protomol/type/Real.h>
%include <protomol/base/Report.h>
%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/bonded/MTorsionSystemForce.h>
%include <protomol/force/bonded/HarmDihedralSystemForce.h>


%extend ProtoMol::HarmDihedralSystemForce<ProtoMol::PeriodicBoundaryConditions> {
   void setPars(Real k, int dih, double dihref) {
      std::vector<Value> v;
      v.push_back(Value(k));
      v.push_back(Value(dih));
      v.push_back(Value(dihref));
      v.push_back(Value(false));
      self->doSetParameters(v);
   }
}

%extend ProtoMol::HarmDihedralSystemForce<ProtoMol::VacuumBoundaryConditions> {
   void setPars(Real k, int dih, double dihref) {
      std::vector<Value> v;
      v.push_back(Value(k));
      v.push_back(Value(dih));
      v.push_back(Value(dihref));
      v.push_back(Value(false));
      self->doSetParameters(v);
   }
}


%template(MTSF_Periodic) ProtoMol::MTorsionSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(MTSF_Vacuum) ProtoMol::MTorsionSystemForce <ProtoMol::VacuumBoundaryConditions >;
%template(HDSF_Periodic) ProtoMol::HarmDihedralSystemForce<ProtoMol::PeriodicBoundaryConditions >;
%template(HDSF_Vacuum) ProtoMol::HarmDihedralSystemForce<ProtoMol::VacuumBoundaryConditions >;
