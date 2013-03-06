%module GenericTopology
%{
//#include "CubicCellManager.h"
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/topology/Topology.h>
#include <protomol/topology/BuildTopology.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>
#include <protomol/topology/LennardJonesParameters.h>
#include <protomol/type/PSF.h>
#include <protomol/type/PAR.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_vector.i"
%include <protomol/type/Vector3DBlock.h>
%include <protomol/topology/Atom.h>
%include <protomol/topology/Bond.h>
%include <protomol/topology/Angle.h>
%include <protomol/topology/Torsion.h>
%include <protomol/topology/GenericTopology.h>
%include <protomol/topology/SemiGenericTopology.h>
%include <protomol/topology/BuildTopology.h>
%include <protomol/topology/LennardJonesParameters.h>
%template(SGT_Periodic) ProtoMol::SemiGenericTopology <ProtoMol::PeriodicBoundaryConditions >;
%template(SGT_Vacuum) ProtoMol::SemiGenericTopology <ProtoMol::VacuumBoundaryConditions >;

%include <protomol/topology/Topology.h>

%template(T_Periodic) ProtoMol::Topology <ProtoMol::PeriodicBoundaryConditions, ProtoMol::CubicCellManager >;
%template(T_Vacuum) ProtoMol::Topology <ProtoMol::VacuumBoundaryConditions, ProtoMol::CubicCellManager >;


%extend ProtoMol::GenericTopology {
        void setExclusion(char* e) {self->exclude = ExclusionType((const char*)e);}
        CoulombSCPISMParameterTable* makeSCPISM() {CoulombSCPISMParameterTable* SCPISMParameters = new CoulombSCPISMParameterTable(); SCPISMParameters->populateTable(); return SCPISMParameters;}
};
%extend ProtoMol::Topology<ProtoMol::VacuumBoundaryConditions, ProtoMol::CubicCellManager> {
	void resizeLennardJonesParameters(int sz) {
           self->lennardJonesParameters.resize(sz);
        }
	void setLennardJonesParameters(int type1, int type2, ProtoMol::LennardJonesParameters params) {
           self->lennardJonesParameters.set(type1, type2, params);
        }
        void setCellSize(ProtoMol::Real r){
	   self->cellManager.setCellSize(r);
	}
        Atom getAtom(int i) {return self->atoms[i];}
        Bond getBond(int i) {return self->bonds[i];}
        Angle getAngle(int i) {return self->angles[i];}
        Torsion getDihedral(int i) {return self->dihedrals[i];}
        Torsion getImproper(int i) {return self->impropers[i];}
        
}

%extend ProtoMol::Topology<ProtoMol::PeriodicBoundaryConditions, ProtoMol::CubicCellManager> {
	void setCellSize(ProtoMol::Real r){
	   self->cellManager.setCellSize(r);
	}
	void resizeLennardJonesParameters(int sz) {
           self->lennardJonesParameters.resize(sz);
        }
	void setLennardJonesParameters(int type1, int type2, ProtoMol::LennardJonesParameters params) {
           self->lennardJonesParameters.set(type1, type2, params);
        }

        void setBC(float cB1x, float cB1y, float cB1z,
                   float cB2x, float cB2y, float cB2z,
                   float cB3x, float cB3y, float cB3z,
                   float cOx, float cOy, float cOz) {
           self->boundaryConditions.set(ProtoMol::Vector3D(cB1x, cB1y, cB1z),
                                        ProtoMol::Vector3D(cB2x, cB2y, cB2z),
                                        ProtoMol::Vector3D(cB3x, cB3y, cB3z),
                                        ProtoMol::Vector3D(cOx, cOy, cOz));
        }
        Atom getAtom(int i) {return self->atoms[i];}
        Bond getBond(int i) {return self->bonds[i];}
        Angle getAngle(int i) {return self->angles[i];}
        Torsion getDihedral(int i) {return self->dihedrals[i];}
        Torsion getImproper(int i) {return self->impropers[i];}
};
		
		
	
