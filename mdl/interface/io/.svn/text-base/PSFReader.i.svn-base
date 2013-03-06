%module PSFReader
%{
#include <protomol/io/PSFReader.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include <protomol/type/PSF.h>
%include <protomol/io/PSFReader.h>


%extend ProtoMol::PSF
{
   PSF::Atom getAtom(int i) {return self->atoms[i];}
   PSF::Bond getBond(int i) {return self->bonds[i];}
   PSF::Angle getAngle(int i) {return self->angles[i];}
   PSF::Dihedral getDihedral(int i) {return self->dihedrals[i];}
   PSF::Improper getImproper(int i) {return self->impropers[i];}
   PSF::Donor getDonor(int i) {return self->donors[i];}
   PSF::Acceptor getAcceptor(int i) {return self->acceptors[i];}
   PSF::Nonbonded getNonbonded(int i) {return self->nonbondeds[i];}
   PSF::Ngrp getNgrp(int i) {return self->ngrp[i];}
   Real getMass(int i) {return self->atoms[i].mass;}

   // A bit ugly right now, but SWIG did not support nested structs
   // and I do not want to change ProtoMol 3 to unnest everything.
   int getAttributeInt(std::string structure, int i, std::string attribute) {
      if (structure == "atom") {
         if (attribute == "number") return self->atoms[i].number;
         else if (attribute == "residue_sequence") return self->atoms[i].residue_sequence;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }
      else if (structure == "bond") {
         if (attribute == "number") return self->bonds[i].number;
         else if (attribute == "atom1") return self->bonds[i].atom1;
         else if (attribute == "atom2") return self->bonds[i].atom2;
	 else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }
      else if (structure == "angle") {
         if (attribute == "number") return self->angles[i].number;
         else if (attribute == "atom1") return self->angles[i].atom1;
         else if (attribute == "atom2") return self->angles[i].atom2;
	 else if (attribute == "atom3") return self->angles[i].atom3;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }
      else if (structure == "dihedral") {
         if (attribute == "number") return self->dihedrals[i].number;
         else if (attribute == "atom1") return self->dihedrals[i].atom1;
         else if (attribute == "atom2") return self->dihedrals[i].atom2;
	 else if (attribute == "atom3") return self->dihedrals[i].atom3;
	 else if (attribute == "atom4") return self->dihedrals[i].atom4;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }      
      else if (structure == "improper") {
         if (attribute == "number") return self->impropers[i].number;
         else if (attribute == "atom1") return self->impropers[i].atom1;
         else if (attribute == "atom2") return self->impropers[i].atom2;
	 else if (attribute == "atom3") return self->impropers[i].atom3;
	 else if (attribute == "atom4") return self->impropers[i].atom4;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }
      else if (structure == "donor") {
         if (attribute == "number") return self->donors[i].number;
         else if (attribute == "atom1") return self->donors[i].atom1;
         else if (attribute == "atom2") return self->donors[i].atom2;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }
      else if (structure == "acceptor") {
         if (attribute == "number") return self->acceptors[i].number;
         else if (attribute == "atom1") return self->acceptors[i].atom1;
         else if (attribute == "atom2") return self->acceptors[i].atom2;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         }
      }   
      else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
      }   
   }

   Real getAttributeReal(std::string structure, int i, std::string attribute)  {
      if (structure == "atom") {
         if (attribute == "charge") return self->atoms[i].charge;
         else if (attribute == "mass") return self->atoms[i].mass;
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
         } 
      }
      else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return -1;
      }   
   }

   char* getAttributeString(std::string structure, int i, std::string attribute)   {
      if (structure == "atom") {
         if (attribute == "seg_id") return (char*)(self->atoms[i].seg_id.c_str());
         else if (attribute == "residue_name") return (char*)(self->atoms[i].residue_name.c_str());
         else if (attribute == "atom_name") return (char*)(self->atoms[i].atom_name.c_str());
         else if (attribute == "atom_type") return (char*)(self->atoms[i].atom_type.c_str());
         else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return "";
         } 
      }
      else {
            cout << "INVALID ATTRIBUTE " << attribute << endl;
            return "";
      }  
   }

   int numAtoms() {return self->atoms.size();}
   int numBonds() {return self->bonds.size();}
   int numAngles() {return self->angles.size();}
   int numDihedrals() {return self->dihedrals.size();}
   int numImpropers() {return self->impropers.size();}
   int numDonors() {return self->donors.size();}
   int numAcceptors() {return self->acceptors.size();}
   int numNonbondeds() {return self->nonbondeds.size();}
   int numNgrps() {return self->ngrp.size();}
}