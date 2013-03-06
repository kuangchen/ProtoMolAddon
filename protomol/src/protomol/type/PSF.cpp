#include <protomol/type/PSF.h>

using namespace ProtoMol;
//____ PSF

void PSF::clear() {
  atoms.clear();
  bonds.clear();
  angles.clear();
  dihedrals.clear();
  impropers.clear();
  donors.clear();
  acceptors.clear();
  nonbondeds.clear();
  ngrp.clear();

  //Amber style dihedrals
  rb_dihedrals.clear();

}

