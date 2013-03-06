#include <protomol/integrator/hessian/Hessian.h>

#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/force/hessian/ReducedHessAngle.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/C1SwitchingFunction.h>

#include <protomol/force/hessian/HessDihedral.h>
#include <protomol/force/hessian/ReducedHessBond.h>
#include <protomol/force/hessian/ReducedHessCoulomb.h>
#include <protomol/force/hessian/ReducedHessCoulombDiElec.h>
#include <protomol/force/hessian/ReducedHessCoulombSCPISM.h>
#include <protomol/force/hessian/ReducedHessLennardJones.h>
#include <protomol/force/hessian/ReducedHessBornSelf.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/born/BornRadii.h>

#include <protomol/force/GB/GBBornBurialTerm.h>
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/force/GB/GBACEForce.h>
#include <protomol/force/hessian/ReducedHessGBACE.h>
//#include <protomol/force/hessian/GB/ReducedHessGBForce.h>

#include <protomol/force/GB/GBForce.h>
#include <protomol/force/hessian/ReducedHessGB.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

namespace ProtoMol {
//constructors
Hessian::Hessian() {
  hessM = 0;
}

Hessian::Hessian(unsigned int szin) {
  sz = szin;
  try{
    hessM = new double[sz * sz];   //assign array
  }catch(bad_alloc&){
    report << error << "[Hessian::Hessian] Cannot allocate memory for Hessian!" << endr;
  }
}

Hessian::~Hessian() {
  if (hessM != 0) delete[] hessM;
}

// copy constructor
Hessian::Hessian(const Hessian &hess) {
  sz = hess.sz;
  if (hess.hessM != 0) {
    try{
        hessM = new double[sz * sz];   //assign array
    }catch(bad_alloc&){
        report << error << "[Hessian::Hessian] Cannot allocate memory for Hessian!" << endr;
    }
    //hessM = new double[sz * sz];
    for (unsigned int i = 0; i < sz * sz; i++) hessM[i] = hess.hessM[i];
  } else
    hessM = 0;
  myBond = hess.myBond;
  myAngle = hess.myAngle;
  myCoulomb = hess.myCoulomb;
  myCoulombDielec = hess.myCoulombDielec;
  myCoulombSCPISM = hess.myCoulombSCPISM;
  myBornRadii = hess.myBornRadii;
  myBornSelf = hess.myBornSelf;
  //GB
  myGBBornBurialTerm = hess.myGBBornBurialTerm;
  myGBBornRadii = hess.myGBBornRadii;
  myGBACEForce = hess.myGBACEForce;
  myGBForce = hess.myGBForce;

  //testing only

  myLennardJones = hess.myLennardJones;
  myDihedral = hess.myDihedral;
  myRBDihedral = hess.myRBDihedral;
  myImproper = hess.myImproper;
  cCutoff = hess.cCutoff;
  cSwitchon = hess.cSwitchon;
  cSwitch = hess.cSwitch;
  cOrder = hess.cOrder;
  cSwitchoff = hess.cSwitchoff;
  lCutoff = hess.lCutoff;
  lSwitchon = hess.lSwitchon;
  lSwitch = hess.lSwitch;
  lOrder = hess.lOrder;
  lSwitchoff = hess.lSwitchoff;
}

void Hessian::initialData(unsigned int szin) {
  sz = szin;
  if (hessM == 0){
      try{
        hessM = new double[sz * sz];   //assign array
      }catch(bad_alloc&){
        report << error << "[Hessian::initialData] Cannot allocate memory for Hessian!" << endr;
      }
  }
  //if (hessM == 0) hessM = new double[sz * sz]; //assign array
}

void Hessian::findForces(ForceGroup *overloadedForces) {
  vector<Force *> ListForces = overloadedForces->getForces();
  //
  lCutoff = lSwitchon = lSwitch = cCutoff = cSwitchon = cSwitch = 0.0;
  lOrder = cOrder = lSwitchoff = cSwitchoff = 0.0;
  Real mCutoff = 0.0;
  D = 78.0; S = 0.3; epsi = 1.0;
  myBornSwitch = 3; myDielecConst = 80.0;
  myBond = 
      myAngle = 
          myCoulomb = 
              myCoulombDielec = 
                   myCoulombSCPISM =
                       myLennardJones =
                           myDihedral =
                              myImproper = 
                                  myBornRadii =
                                       myBornSelf = 
                                           myRBDihedral = 
                                               myGBBornBurialTerm = 
                                                   myGBBornRadii = 
                                                       myGBACEForce = 
                                                                myGBForce = false;

   solvationparam = 2.26 / 418.4 ; watersphereradius = 1.4;

  for (unsigned int i = 0; i < ListForces.size(); i++) {

    if (equalNocase(ListForces[i]->getId(), "Bond")) {
      myBond = true;
    } else if (equalNocase(ListForces[i]->getId(), "Angle")) {
      myAngle = true;
    // Check for combined force
    } else if (equalStartNocase("LennardJones CoulombSCPISM BornRadii", ListForces[i]->getId())) {
      myBornRadii = true; myLennardJones = true; myCoulombSCPISM = true;
      vector<Parameter> Fparam;
      lSwitch = cSwitch = 1;
      int curr_force = 0;
      bool last_was_cutoff = false;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          //get outer cuttoff for C1/C2 switches
          if(curr_force == 0){
            cCutoff = Fparam[j].value;
          }else if(curr_force == 1){          
            lCutoff = Fparam[j].value;
          }
          if(!last_was_cutoff){
            curr_force++;
          }else{
            //if two cuttoffs in sequence then outer value
            mCutoff = Fparam[j].value;
          }
          last_was_cutoff = true;
        } else {
          last_was_cutoff = false;
          if (equalNocase(Fparam[j].keyword, "-switchon")) {
            if(curr_force == 0){
              cSwitchon = Fparam[j].value;
              if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
              else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
              else cSwitch = 1;
            }else if(curr_force == 1){
              lSwitchon = Fparam[j].value;
              if (equalStartNocase("C2", Fparam[j].text)) lSwitch = 2;
              else if (equalStartNocase("Cn", Fparam[j].text)) lSwitch = 3;
              else lSwitch = 1;
            }
          } else if (equalNocase(Fparam[j].keyword, "-n")) {
            if(curr_force == 0){
              cOrder = Fparam[j].value;
              cSwitch = 3;
            }else if(curr_force == 1){          
              lOrder = Fparam[j].value;
              lSwitch = 3;
            }
          } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
            if(curr_force == 0){
              cSwitchoff = Fparam[j].value;
              cSwitch = 3;
            }else if(curr_force == 1){          
              lSwitchoff = Fparam[j].value;
              lSwitch = 3;
            }
          }
        }
      }
    // Single forces here
    } else if (equalStartNocase("CoulombDiElec", ListForces[i]->getId())) {
      myCoulombDielec = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text))
            cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-D"))
          D = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-S"))
          S = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-EPS"))
          epsi = Fparam[j].value;
      }
    } else if (equalStartNocase("CoulombSCPISM", ListForces[i]->getId())) {
      myCoulombSCPISM = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("Coulomb", ListForces[i]->getId())) {
      myCoulomb = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("LennardJones", ListForces[i]->getId())) {
      myLennardJones = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          lCutoff = Fparam[j].value;
          if (lSwitch == 0) lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          lSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) lSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) lSwitch = 3;
          else lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          lOrder = Fparam[j].value;
          lSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          lSwitchoff = Fparam[j].value;
          lSwitch = 3;
        }
      }
    } else if (equalStartNocase("Dihedral", ListForces[i]->getId())) {
      myDihedral = true;
    } else if (equalStartNocase("RBDihedral", ListForces[i]->getId())) {
      myRBDihedral = true;
    } else if (equalStartNocase("Improper", ListForces[i]->getId())) {
      myImproper = true;
    } else if (equalStartNocase("BornRadii", ListForces[i]->getId())) {
      myBornRadii = true;
    } else if (equalStartNocase("BornSelf", ListForces[i]->getId())) {
      myBornSelf = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-bornswitch")) {
          myBornSwitch = Fparam[j].value;
        } else if (equalNocase(Fparam[j].keyword, "-D")) {
          myDielecConst = Fparam[j].value;
        } 
      }
    } else if (equalStartNocase("GBBornBurialTerm", ListForces[i]->getId())) {
        myGBBornBurialTerm = true;
    } else if (equalStartNocase("GBBornRadii", ListForces[i]->getId())) {
        myGBBornRadii = true;
    } else if (equalStartNocase("GBACEForce", ListForces[i]->getId())) {
        myGBACEForce = true;
        vector<Parameter> Fparam;
        ListForces[i]->getParameters(Fparam);
        for (unsigned int j = 0; j < Fparam.size(); j++) {
           if (equalNocase(Fparam[j].keyword, "-solvationparam")) {
              solvationparam = Fparam[j].value;
           }else if (equalNocase(Fparam[j].keyword,"-watersphereradius")) {
              watersphereradius = Fparam[j].value;
           }
        }
    } else if (equalStartNocase("GBForce", ListForces[i]->getId())) {

      myGBForce = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
         if (equalNocase(Fparam[j].keyword, "-soluteDielec")) {
            soluteDielec = Fparam[j].value;
         }else if (equalNocase(Fparam[j].keyword,"-solventDielec")) {
            solventDielec = Fparam[j].value;
         }
      }
    }

  }      
  //set maxiumum cutoff value
  Real tempCoff = max(cCutoff, lCutoff);
  //allow for 'outer' cutoff
  tempCoff = max(mCutoff, tempCoff);

  cutOff = max(cSwitchoff, lSwitchoff);
  cutOff = max(cutOff, tempCoff);
}

void Hessian::evaluate(const Vector3DBlock *myPositions,
                       GenericTopology *myTopo,
                       bool mrw) {
  int a1, a2, a3;
  unsigned int i;
  ReducedHessAngle rh;
  Matrix3By3 rha;

  sz = 3 * myPositions->size();
  //Impropers
  if (myImproper) {
    HessDihedral hi;        //create improper hessian
    for (unsigned int i = 0; i < myTopo->impropers.size(); i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->impropers[i].multiplicity; j++)
        if (myTopo->impropers[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = myTopo->impropers[i].atom1; aout[1] =
          myTopo->impropers[i].atom2;
        aout[2] = myTopo->impropers[i].atom3; aout[3] =
          myTopo->impropers[i].atom4;
        //
        hi.evaluate(myTopo->impropers[i], (*myPositions)[aout[0]],
                    (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                    (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++)
          for (int kk = 0; kk < 4; kk++) {
            Matrix3By3 rhd = hi(ii, kk);
            outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
                myTopo->atoms[aout[kk]].scaledMass, rhd, mrw, sz, hessM);
          }

      }
    }
  }

  //Dihedrals
  if (myDihedral) {
    HessDihedral hd;        //create dihedral hessian
    for (unsigned int i = 0; i < myTopo->dihedrals.size(); i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->dihedrals[i].multiplicity; j++)
        if (myTopo->dihedrals[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = myTopo->dihedrals[i].atom1; aout[1] =
          myTopo->dihedrals[i].atom2;
        aout[2] = myTopo->dihedrals[i].atom3; aout[3] =
          myTopo->dihedrals[i].atom4;
        //
        hd.evaluate(myTopo->dihedrals[i], (*myPositions)[aout[0]],
                    (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                    (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++)
          for (int kk = 0; kk < 4; kk++) {
            Matrix3By3 rhd = hd(ii, kk);
            outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
                myTopo->atoms[aout[kk]].scaledMass, rhd, mrw, sz, hessM);
          }

      }
    }
  }

  //RBDihedrals
  if (myRBDihedral) {
    HessDihedral hd;        //create dihedral hessian
    for (unsigned int i = 0; i < myTopo->rb_dihedrals.size(); i++) {
      RBTorsion rbt = myTopo->rb_dihedrals[i];
      bool nonZForce = false;       //test for force constants

      if (rbt.C0 || rbt.C1 || rbt.C2 || rbt.C3 || rbt.C4 || rbt.C5) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = rbt.atom1; aout[1] = rbt.atom2;
        aout[2] = rbt.atom3; aout[3] = rbt.atom4;
        //
        hd.evaluate(rbt, (*myPositions)[aout[0]],
                    (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                    (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++)
          for (int kk = 0; kk < 4; kk++) {
            Matrix3By3 rhd = hd(ii, kk);
            outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
                myTopo->atoms[aout[kk]].scaledMass, rhd, mrw, sz, hessM);
          }

      }

    }
  }

  //Bonds
  if (myBond){
    for (i = 0; i < myTopo->bonds.size(); i++) {
      a1 = myTopo->bonds[i].atom1; a2 = myTopo->bonds[i].atom2;
      Real r_0 = myTopo->bonds[i].restLength;
      Real k = myTopo->bonds[i].springConstant;
      Matrix3By3 bondHess12 =
        reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
      //output sparse matrix
      outputSparsePairMatrix(a1,a2,myTopo->atoms[a1].scaledMass,myTopo->atoms[a2].scaledMass,
                                bondHess12,mrw,sz,hessM);

    }
  }

  //Angles
    if (myAngle){
    for (i = 0; i < myTopo->angles.size(); i++) {
      a1 = myTopo->angles[i].atom1;
      a2 = myTopo->angles[i].atom2;
      a3 = myTopo->angles[i].atom3;
      Real theta0 = myTopo->angles[i].restAngle;
      Real k_t = myTopo->angles[i].forceConstant;
      Real ubConst = myTopo->angles[i].ureyBradleyConstant;
      Real ubRestL = myTopo->angles[i].ureyBradleyRestLength;
      // ReducedHessAngle for atoms a1, a2 and a3
      rh.evaluate((*myPositions)[a1], (*myPositions)[a2], (*myPositions)[a3],
                  k_t,
                  theta0);
      //ureyBradley
      if (ubConst) {
        //Cheat using bond hessian as same as UB!!!!
        Matrix3By3 ubm =
          reducedHessBond((*myPositions)[a1], (*myPositions)[a3], ubConst,
                          ubRestL);
        rh.accumulateTo(0, 0, ubm);
        rh.accumulateTo(2, 2, ubm);
        rh.accumulateNegTo(2, 0, ubm);
        rh.accumulateNegTo(0, 2, ubm);
      }
      //output sparse matrix
      int aout[3];
      aout[0] = a1; aout[1] = a2; aout[2] = a3;
      for (int ii = 0; ii < 3; ii++)
        for (int kk = 0; kk < 3; kk++) {
          rha = rh(ii, kk);
          outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
              myTopo->atoms[aout[kk]].scaledMass, rha, mrw, sz, hessM);

        }

    }
  }
  //Pairwise forces

  //Pre-calculate Born radii if Self energy Hessian required
  if(myBornRadii && myBornSelf && myTopo->doSCPISM){
    evaluateBornRadii(myPositions, myTopo);
  }

  if (myGBBornBurialTerm && myTopo->doGBSAOpenMM) {
     //report << plain <<"Hessian : Calculate forces GBBornBurialTerm"<<endr;
     evaluateGBBornBurialTerm(myPositions, myTopo);

/*
     for(int i=0; i<myTopo->atoms.size();i++)
        for (int j = i+1 ; j<myTopo->atoms.size() ; j++) 
            report << plain <<"From force : Atom1 "<<i<<", Atom2 "<<j<<" dist "<<myTopo->atoms[i].myGBSA_T->distij[j]<<" Lij "<<myTopo->atoms[i].myGBSA_T->Lvalues[j]<<endr;
*/
  }

  if (myGBBornRadii && myTopo->doGBSAOpenMM) {
     //report << plain <<"Hessian : Calculate forces GBBornRadii"<<endr;
     evaluateGBBornRadii(myPositions, myTopo);
  }

  if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Appropriate flags set for calculation of GB hessian"<<endr;
  }

  unsigned int atoms_size = myTopo->atoms.size();
  for (unsigned int i = 0; i < atoms_size; i++){
    for (unsigned int j = i + 1; j < atoms_size; j++){ 
      Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      //Lennard jones
      if (myLennardJones) 
        rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, mrw);
      //Coulombic
      if (myCoulomb)
        rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, mrw);
      //Coulombic Implicit solvent
      if (myCoulombDielec)
        rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, mrw);
      //SCP
      if (myCoulombSCPISM)
        rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, mrw);
      //Bourn radii
      if (myBornRadii && myBornSelf && myTopo->doSCPISM)          
        rhp += evaluateBornSelfPair(i, j, myPositions, myTopo);


      if (myGBBornBurialTerm && myGBBornRadii && myGBACEForce && myTopo->doGBSAOpenMM) {
        rhp += evaluateGBACEPair(i, j, myPositions, myTopo);
      }

      if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
         rhp += evaluateGBPair(i, j, myPositions, myTopo);
      }
      //output sum to matrix
      outputSparsePairMatrix(i, j, myTopo->atoms[i].scaledMass, myTopo->atoms[j].scaledMass,
                              rhp, mrw, sz, hessM);
      //
    }
  }
  //

}

void Hessian::evaluatePairs(int i, int j, int pairType, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo, bool mrw, 
                                            int mat_i, int mat_j, int mat_sz, double * mat_array) {

  Matrix3By3 rha = evaluatePairsMatrix(i, j, pairType, myPositions, myTopo, mrw);
  //
  outputSparsePairMatrix(mat_i,mat_j,myTopo->atoms[i].scaledMass,myTopo->atoms[j].scaledMass,
                          rha,mrw,mat_sz,mat_array);

}

Matrix3By3 Hessian::evaluatePairsMatrix(int i, int j, int pairType, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo, bool mrw) {

  Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  ExclusionClass ec;
  //Coulomb
  ReducedHessCoulomb rHess;
  CoulombForce hForce;
  ReducedHessCoulombDielec rHessDi;
  CoulombForceDiElec hForceDi;
  ReducedHessCoulombSCPISM rHessScp;
  CoulombSCPISMForce hForceScp;
  //LJ
  ReducedHessLennardJones rHessLj;
  LennardJonesForce hForceLj;
  //SCP
  Real alpha_ij = 0;
  //Switches
  Real pCutoff, pSwitchon, pSwitch, pOrder, pSwitchoff; 


  if(pairType == LENNARDJONES){	//setup switch paramiters
    pCutoff = lCutoff; pSwitchon = lSwitchon; 
    pSwitch = lSwitch; pOrder = lOrder; pSwitchoff = lSwitchoff; 
  }else{
    pCutoff = cCutoff; pSwitchon = cSwitchon; 
    pSwitch = cSwitch; pOrder = cOrder; pSwitchoff = cSwitchoff; 
  }
  //if not bonded/dihedral
  ec = myTopo->exclusions.check(i, j);
  if (ec != EXCLUSION_FULL) {
      Vector3D rij =
            myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
      //####FAST BAIL ON CUTOFF CHECK?#################################################
      Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      Real a = rij.normSquared();
      Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
      //SCP extras
      if(pairType == COULOMBSCPISM){
          if (!myTopo->atoms[i].mySCPISM_A || !myTopo->atoms[j].mySCPISM_A)
                    report << error << "[Hessian::evaluateCoulombSCPISM] SCPISM data not set." << endr;
          alpha_ij = myTopo->atoms[i].mySCPISM_A->sqrtalphaSCPISM *
                            myTopo->atoms[j].mySCPISM_A->sqrtalphaSCPISM;
      }
      //
      if (pSwitch) {
          if (pSwitch == 3) {
              CnSwitchingFunction cnsf(pOrder, pSwitchon, pSwitchoff, pCutoff);
              cnsf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                  }
                  mz = cnsf.hessian(rij, a);
              }
          } else if (pSwitch == 2) {
              C2SwitchingFunction c2sf(pSwitchon, pCutoff);
              c2sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                break;
                  }
                  mz = c2sf.hessian(rij, a);
              }
          } else {
              C1SwitchingFunction c1sf(pCutoff);
              c1sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                break;
                  }
                  mz = c1sf.hessian(rij, a);
              }
          }
      }
      switch(pairType){
          case LENNARDJONES:  rha = rHessLj(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec);
                              break;
          case COULOMB:		rha = rHess(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec);
                              break;
          case COULOMBDIELEC:	rha = rHessDi(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec, D, S, epsi);
                              break;
          case COULOMBSCPISM: rha = rHessScp(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec, alpha_ij, 80.0);
                              break;
      }
      //
  }
  return rha;
  //
}

void Hessian::outputSparsePairMatrix(int i, int j, Real massi, Real massj, 
                                        Matrix3By3 rha, bool mrw, int arrSz, double *basePoint){
    int eye, jay;
    Real tempf, ms1, ms2, ms3;

    //output sparse matrix
    for (int ll = 0; ll < 3; ll++){
        for (int mm = 0; mm < 3; mm++) {
            eye = i * 3 + 1; jay = j * 3 + 1;
            tempf = rha(ll, mm);
            if (mrw) {
                ms1 = tempf / sqrt(massi * massi);
                ms2 = tempf / sqrt(massj * massj);
                ms3 = -tempf / sqrt(massi * massj);
            }else{
                ms2 = ms1 = tempf;
                ms3 = -tempf;
            }
            basePoint[(eye + ll - 1) + (eye + mm - 1) * arrSz] += ms1;
            basePoint[(jay + ll - 1) + (jay + mm - 1) * arrSz] += ms2;
            basePoint[(eye + ll - 1) + (jay + mm - 1) * arrSz] += ms3;
            basePoint[(jay + ll - 1) + (eye + mm - 1) * arrSz] += ms3;
        }
    }
}

void Hessian::outputSparseMatrix(int i, int j, Real massi, Real massj, 
                                    Matrix3By3 rha, bool mrw, int arrSz, double *basePoint){

    for (int ll = 0; ll < 3; ll++)
        for (int mm = 0; mm < 3; mm++)
            if (mrw)
                basePoint[(i * 3 + ll) + (j * 3 + mm) * arrSz] +=
                        ((rha(ll, mm) / sqrt(massi * massj)));
            else 
                basePoint[(i * 3 + ll) + (j * 3 + mm) * arrSz] +=
                        rha(ll, mm);

}

//Find Born radii
void Hessian::evaluateBornRadii(const Vector3DBlock *myPositions,
                                  GenericTopology *myTopo){

  unsigned int atoms_size = myTopo->atoms.size();
  BornRadii br;
  Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  //SCPISM pre-calculate initialize
  for(unsigned int i=0;i<atoms_size;i++){
    myTopo->atoms[i].mySCPISM_A->bornRadius = myTopo->atoms[i].mySCPISM_A->zeta;
    myTopo->atoms[i].mySCPISM_A->energySum = true;
    myTopo->atoms[i].mySCPISM_A->D_s = 0.0;
  }
  //check all pairs
  for (unsigned int i = 0; i < atoms_size; i++){
    for (unsigned int j = i + 1; j < atoms_size; j++){ 
      //if not bonded/dihedral
      ExclusionClass ec = myTopo->exclusions.check(i, j);
      if (ec != EXCLUSION_FULL) {
        Vector3D rij =
              myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
        Real a = rij.normSquared();
        if(a < BORNCUTOFF2){ //within cutoff of 5A?      
          Real rawE = 0.0, rawF = 0.0;
          br(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);  //do the calculation, no force/energy returned.
        }
      }
    }
  }
}
  
//SCPISM Pairwise Born Self energy Hessian
Matrix3By3 Hessian::evaluateBornSelfPair(int i, int j, const Vector3DBlock *myPositions,
                                       const GenericTopology *myTopo) {

  ReducedHessBornSelf rHessBS;  
  Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  //find the Hessian
  //if not bonded/dihedral
  ExclusionClass ec = myTopo->exclusions.check(i, j);
  if (ec != EXCLUSION_FULL) {
    Vector3D rij =
          myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
    Real a = rij.normSquared();
    if(a < BORNCUTOFF2){ //within cutoff of 5A?      
      rha = rHessBS(a, rij, myTopo, i, j, myBornSwitch, myDielecConst, ec);  //find Hessian
    }
  }
  return rha;
}

//GB Burial term
void Hessian::evaluateGBBornBurialTerm(const Vector3DBlock *myPositions,
                                         GenericTopology *myTopo) {

   unsigned int atom_size = myTopo->atoms.size();
   GBBornBurialTerm gbBornBurialTerm;
   
   //substitute code for preForce() 
   for (unsigned int i = 0 ; i <atom_size ; i++) {
      myTopo->atoms[i].myGBSA_T->preForce();
   }
   for (unsigned int i=0; i < atom_size ; i++) {
      for (unsigned int j = i + 1; j < atom_size; j++){
         ExclusionClass ec = myTopo->exclusions.check(i, j);
         //Now include ALL atoms, no exclusions
         if (1) { //ec != EXCLUSION_FULL) {
           Vector3D rij =
              myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
           Real a = rij.normSquared();
           Real rawE, rawF;
           gbBornBurialTerm(rawE, rawF, a, 1.0/a, rij, myTopo, i, j, ec);
         }
      }
   }
}

//GB born radii
void Hessian::evaluateGBBornRadii(const Vector3DBlock *myPositions,
                                         GenericTopology *myTopo) {

   unsigned int atom_size = myTopo->atoms.size();
   GBBornRadii gbBornRadii;

   for (unsigned int i=0; i < atom_size ; i++) {
      for (unsigned int j = i + 1; j < atom_size; j++){
         ExclusionClass ec = myTopo->exclusions.check(i, j);
         //Now include ALL atoms, no exclusions
         if (1) { //ec != EXCLUSION_FULL) {
           Vector3D rij =
              myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
           Real a = rij.normSquared();
           Real rawE, rawF;
           gbBornRadii(rawE, rawF, a, 1.0/a, rij, myTopo, i, j, ec);
         }
      }
   }
}

//GB ACE force
Matrix3By3 Hessian::evaluateGBACEPair(int i, int j, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo) {

   ReducedHessGBACE rHessGBACE;
   Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

   ExclusionClass ec = myTopo->exclusions.check(i, j);
   //Now include ALL atoms, no exclusions
   if (1) { //ec != EXCLUSION_FULL) {
      Vector3D rij = myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
      Real a = rij.normSquared();
      rha = rHessGBACE(a, rij, myTopo, i, j, solvationparam, watersphereradius, ec);
   }

   return rha;
}

//GB force
Matrix3By3 Hessian::evaluateGBPair(int i, int j, const Vector3DBlock *myPositions, 
                                  const GenericTopology *myTopo) {

   ReducedHessGB rHessGB;
   Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

   int sz = (int)myTopo->atoms.size();
   ExclusionClass ec = myTopo->exclusions.check(i, j);
   //Now include ALL atoms, no exclusions
   if (1) { //ec != EXCLUSION_FULL) {
      Vector3D rij = myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
      Real a = rij.normSquared();
      rha = rHessGB(a, rij, myTopo, i, j, sz, soluteDielec, solventDielec, ec);
   }

   return rha;

}

//Clear Hessian
void Hessian::clear() {
  if (hessM != 0)
    for (unsigned int i = 0; i < sz * sz; i++) hessM[i] = 0.0;

}
    //set Hessian column
    bool Hessian::setHessianColumn( const Vector3DBlock &hescol, const unsigned int columnNumber, 
                                    const GenericTopology *myTopo, const bool massWeight ){
        //get size
        sz = 3 * hescol.size();

        //loop over column elements
        if (hessM != 0 && columnNumber < sz){
            for (unsigned int i = 0; i < sz; i++){
                
                Real factor = 1.0;
                
                if(massWeight) factor /= sqrt(myTopo->atoms[i/3].scaledMass);
                
                hessM[columnNumber * sz + i] += hescol[i/3][i%3] * factor;
            }
            
            return true;
            
        }else{
            return false;
        }
        
    }
    
}
