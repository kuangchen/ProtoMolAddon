#include <protomol/integrator/hessian/BlockHessian.h>

#include <protomol/base/Report.h>
#include <protomol/force/hessian/HessDihedral.h>
#include <protomol/force/hessian/ReducedHessBond.h>
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
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>

//GB
#include <protomol/force/GB/GBBornBurialTerm.h>
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/force/GB/GBACEForce.h>
#include <protomol/force/hessian/ReducedHessGBACE.h>
#include <protomol/force/GB/GBForce.h>
#include <protomol/force/hessian/ReducedHessGB.h>

#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/type/BlockMatrix.h>

#include <protomol/force/bonded/BondSystemForce.h>
#include <protomol/force/bonded/AngleSystemForce.h>
#include <protomol/force/bonded/ImproperSystemForce.h>
#include <protomol/force/bonded/DihedralSystemForce.h>
#include <protomol/force/bonded/RBDihedralSystemForce.h>

#include <iostream>
#include <stdio.h>
#include <fstream>

//defines for including SCPISM and GB and
//adding paiwise forces to evaluateResidues()
#define ADDSCPISM
#define ADDGB
//#define BLOCKPAIRWISEINTERACTION
//#define BLOCKSCPISM
//#define BLOCKGB

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//constructors
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BlockHessian::BlockHessian() {
  residues = 0;
  residues_max = 0;
  residues_alpha_c = residues_phi_n = residues_psi_c = 0;
  hess_array_point = 0;
  hess_eig_point = 0;
  atom_residue = 0;
  atom_res_num = 0;
  blocks_max = 0; 
  atom_block = 0;
  atom_block_num = 0;
  sqrtMass = 0;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BlockHessian::~BlockHessian() {		

  if(residues!=0) delete [] residues;
  if(residues_max!=0) delete [] residues_max;
  if(residues_alpha_c!=0) delete [] residues_alpha_c;
  if(residues_phi_n!=0) delete [] residues_phi_n;
  if(residues_psi_c!=0) delete [] residues_psi_c;
  if(hess_array_point!=0) delete [] hess_array_point;
  if(hess_eig_point!=0) delete [] hess_eig_point;
  if(atom_residue!=0) delete [] atom_residue;
  if(atom_res_num!=0) delete [] atom_res_num;	
  if(blocks_max!=0) delete [] blocks_max;
  if(atom_block!=0) delete [] atom_block;
  if(atom_block_num!=0) delete [] atom_block_num;
  if(sqrtMass!=0) delete [] sqrtMass;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::initialResidueData(const GenericTopology *myTopo, int res_per_block, bool fullElect) {
  //Residues
  int _N = myTopo->atoms.size();
  int _3N = 3* _N;
  //set full electrostatic flag
  fullElectrostatics = fullElect;
  //find number residues
  vector<int> residue_id;
  int current_id = myTopo->atoms[0].residue_seq;
  residue_id.push_back(current_id);
  for(int j=0;j<_N;j++){
    if(myTopo->atoms[j].residue_seq != current_id){
      current_id = myTopo->atoms[j].residue_seq;
      residue_id.push_back(current_id);
    }
  }
  num_residues = residue_id.size();
  //Check
  if(res_per_block < 1 || res_per_block >= num_residues)
      report << error << "[BlockHessian::initialResidueData] Residues per block must be greater than 0 and less than " << num_residues << "." << endr;
  //number of blocks
  num_blocks = (int)ceil((float)num_residues / (float)res_per_block);
  rpb = res_per_block;
  blocks_max = new int[num_blocks];        
  //##report << hint << "blocks "<<num_blocks<<", rpb "<<rpb<<" residues"<<num_residues<<endl;
  //assignments
  try{
    //assign residue array
    residues = new int[num_residues * MAX_ATOMS_PER_RES];   
    residues_max = new int[num_residues];
    residues_alpha_c = new int[num_residues];
    residues_phi_n = new int[num_residues];
    residues_psi_c = new int[num_residues];
    //
    hess_array_point = new int[num_blocks];
    hess_eig_point = new int[num_blocks];
    //reverse lookup for atoms
    atom_residue = new int[myTopo->atoms.size()];
    atom_res_num = new int[myTopo->atoms.size()];
    atom_block = new int[myTopo->atoms.size()];
    atom_block_num = new int[myTopo->atoms.size()];
    //pre calculated square roots
    sqrtMass = new Real[myTopo->atoms.size()];
    //
  }catch(bad_alloc&){
      report << error << "[BlockHessian::initialData] Cannot allocate memory for residue variables!" << endr;
  }
  //
  for(int i=0;i<num_residues * MAX_ATOMS_PER_RES;i++) residues[i] = -1; //flag no atom as default
  //find atoms and alpha carbons belonging to residues
  int last_atom = 0;
  for(int i=0;i<num_residues;i++){
      residues_alpha_c[i] = residues_phi_n[i] = residues_psi_c[i] = -1;	//flag invalid (and phi/psi links)
      int res_idx = 0;
      bool in_block = false;
      for(int j=last_atom;j<_N;j++){
          if(myTopo->atoms[j].residue_seq == residue_id[i]){
              in_block = true;  //only allow contigous atoms in 1 residue, numbers may be re-used.
              residues[i*MAX_ATOMS_PER_RES+res_idx] = j;
              if(myTopo->atoms[j].name.compare("CA")==0){	//find alpha carbons
                  residues_alpha_c[i] = j;	//flag valid and mark alpha carbons
              }
              //
              atom_residue[j] = i;
              atom_res_num[j] = res_idx;
              //
              res_idx++;
          }else{
            last_atom = j;
            if(in_block) break; //out if was in block
          }
      }
      residues_max[i] = res_idx;	//set max number of atoms per residue
      if(residues_alpha_c[i] == -1) report << error << "Alpha Carbon not found for residue "<<
          myTopo->atoms[i].residue_name<<", "<<residue_id[i]<<"."<<endr;
      report << debug(3) << "residue "<<residue_id[i]<<", valid "<<residues_alpha_c[i]<<", max "<<residues_max[i]<<
          ", total "<<num_residues<<", alpha "<<myTopo->atoms[residues_alpha_c[i]].name<<endr;
  }
  //Find N,H for \phis and C,O for \psis
  int bonds_size = myTopo->bonds.size();
  for(int i=0;i<bonds_size;i++){	//N and C
      int atom1 = myTopo->bonds[i].atom1;
      int atom2 = myTopo->bonds[i].atom2;
      if(myTopo->atoms[atom1].name.compare("CA")==0){				//atom 1 alpha_c
              if(myTopo->atoms[atom2].name.compare("N")==0)			//and 2 nitrogen then \phi
                  residues_phi_n[atom_residue[atom1]] = atom2;
              if(myTopo->atoms[atom2].name.compare("C")==0)			//and 2 carbon then \psi
                  residues_psi_c[atom_residue[atom1]] = atom2;
      }else{																			//OR
          if(myTopo->atoms[atom2].name.compare("CA")==0){			//atom 2 alpha_c
                  if(myTopo->atoms[atom1].name.compare("N")==0)		//and 1 nitrogen then \phi
                      residues_phi_n[atom_residue[atom2]] = myTopo->bonds[i].atom1;
                  if(myTopo->atoms[atom1].name.compare("C")==0)		//and 1 carbon then \psi
                      residues_psi_c[atom_residue[atom2]] = atom1;
          }
      }
  }
  for(int i=0;i<num_residues;i++){	//diags
      report << debug(3) << "residue "<<i<<", phi N "<<residues_phi_n[i]<<
          ", psi C "<<residues_psi_c[i]<<", alpha C "<<residues_alpha_c[i]<<endr;
      if(residues_phi_n[i] == -1 || residues_psi_c[i] == -1)
          report << error << "Phi/Psi parameter not found!" << endr;
  }
  //test atoms are in the correct order (if not we need to implement eigenvector row swapping)
  int i_res = 0;
  int i_res_idx = 0;
  for(int i=0;i<_N;i++){
      if(residues[i_res*MAX_ATOMS_PER_RES+i_res_idx] != i) 
          report << error << "Atoms NOT in sequence at residue "<<i_res<<", atom "<<i_res_idx
                  <<", number "<<residues[i_res*MAX_ATOMS_PER_RES+i_res_idx]<<", i "<<i<<endr;
      i_res_idx++;
      if(i_res_idx >= residues_max[i_res]){
          i_res++;
          i_res_idx = 0;
      }
  }
  report << hint << "Atoms ARE in sequence! Number of residues " << num_residues << "." << endr;
  //find max of max atoms per residue
  hess_array_size = 0;
  hess_eig_size = 0;
  for(int i=0;i<num_blocks;i++){
      hess_array_point[i] = hess_array_size;
      hess_eig_point[i] = hess_eig_size;
      int block_max = 0;
      for(int j=0;j<res_per_block;j++)
          if(i*res_per_block+j<num_residues)
              block_max += residues_max[i*res_per_block+j];
      blocks_max[i] = block_max;    
      //##report << hint << "Block "<<i<<" max " << block_max << endr;
      hess_array_size += block_max * block_max * 9; 
      hess_eig_size += block_max;
  }
  //find atom blocks and indexes
  for(int j=0;j<_N;j++){
      atom_block[j] = atom_residue[j] / rpb;
      atom_block_num[j] = atom_res_num[j];
      for(int k=0;k<atom_residue[j] % rpb;k++) 
          atom_block_num[j] += residues_max[atom_block[j] * rpb + k];    
  }
  //Assign array
  rsz = hess_array_size;//+(num_blocks*num_blocks * 9);
  //report << hint << "Reduced Hessian size "<<hess_array_size<<endr;
  //Assign Hessian Blocks
  blocks.resize(num_blocks);
  adj_blocks.resize(num_blocks-1);
  memory_base = 0;
  for(int i=0;i<num_blocks;i++){
    int blk_m3_i = blocks_max[i]*3;
    blocks[i].initialize(hess_eig_point[i]*3,hess_eig_point[i]*3,blk_m3_i,blk_m3_i);  //clear block
    blocks[i].clear();  //clear block
    memory_base += blk_m3_i * blk_m3_i;
    if(i<num_blocks - 1){
      int blk_m3_ip1 = blocks_max[i+1]*3;
      adj_blocks[i].initialize(hess_eig_point[i]*3,hess_eig_point[i+1]*3,blk_m3_i,blk_m3_ip1);  //clear block
      adj_blocks[i].clear();  //clear block
      memory_base += blk_m3_i * blk_m3_ip1;
    }
  }
  //pre calculate mass square roots
  for(int j=0;j<_N;j++) sqrtMass[j] = sqrt(myTopo->atoms[j].scaledMass);//
  //full Hessian for electrostatics?
  if(fullElectrostatics){  
    electroStatics.initialize(0,0,_3N,_3N);
    memory_base += _3N * _3N;
  }

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Evaluate course high frequency Hessians
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::evaluateResidues(const Vector3DBlock *myPositions,
                                      GenericTopology *myTopo,
                                        bool simuMin) {
  int a1, a2, a3;
  unsigned int i;
  ReducedHessAngle rh;
  Matrix3By3 rha;
  
  //
  for(int i=0;i<num_blocks;i++)
    blocks[i].clear();  //clear blocks
  
  //RBDihedrals
  if (myRBDihedral){
    
    HessDihedral hd;        //create dihedral hessian
    //loop over all
    for (unsigned int i = 0; i < myTopo->rb_dihedrals.size(); i++) {
      RBTorsion rbt = myTopo->rb_dihedrals[i];
      bool nonZForce = false;       //test for force constants
      
      if (rbt.C0 || rbt.C1 || rbt.C2 || rbt.C3 || rbt.C4 || rbt.C5) nonZForce = true;
      
      if (nonZForce) {
        int aout[4];
        aout[0] = rbt.atom1; aout[1] = rbt.atom2;
        aout[2] = rbt.atom3; aout[3] = rbt.atom4;
        
        //test all in same block
        int ar0 = atom_block[aout[0]];
        if(atom_block[aout[1]] == ar0 && atom_block[aout[2]] == ar0 && atom_block[aout[3]] == ar0){
          //pseudo minimum?
          if(simuMin){
            const Real phi = dihedralAngle(aout, *myPositions);
            //
            RBTorsion currRBTorsion = myTopo->rb_dihedrals[i];
            
            currRBTorsion.Offset = phi;
            

            //
            hd.evaluate(currRBTorsion, (*myPositions)[aout[0]],
                        (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                        (*myPositions)[aout[3]]);
            
          //else normal
          }else{
            //
            hd.evaluate(rbt, (*myPositions)[aout[0]],
                        (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                        (*myPositions)[aout[3]]);
          }
          
          //output sparse matrix
          for (int ii = 0; ii < 4; ii++){
            for (int kk = 0; kk < 4; kk++) {
              Matrix3By3 rhd = hd(ii, kk);
              outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]],
                             sqrtMass[aout[kk]], rhd, blocks[ar0].Rows, blocks[ar0].arrayPointer());
            
            }
          }
                      
        }
        
      }
      
    }
    
  }
  
  //Impropers
  if (myImproper) {
    HessDihedral hi;        //create improper hessian    
    unsigned int impropers_size = myTopo->impropers.size();
    for (unsigned int i = 0; i < impropers_size; i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->impropers[i].multiplicity; j++)
        if (myTopo->impropers[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4] = {myTopo->impropers[i].atom1, myTopo->impropers[i].atom2,
                          myTopo->impropers[i].atom3, myTopo->impropers[i].atom4};
        //test all in same block
        int ar0 = atom_block[aout[0]];
        if(atom_block[aout[1]] == ar0 && atom_block[aout[2]] == ar0 && atom_block[aout[3]] == ar0){
          if(simuMin){
              //    
              Real phi = dihedralAngle(aout, *myPositions);
              //
              Torsion currTorsion = myTopo->impropers[i];
              for (int j = 0; j < currTorsion.multiplicity; j++){
                  currTorsion.phaseShift[j] = phi;
              }
              hi.evaluate(currTorsion, (*myPositions)[aout[0]],
                          (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                          (*myPositions)[aout[3]]);

          }else{
              //
              hi.evaluate(myTopo->impropers[i], (*myPositions)[aout[0]],
                          (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                          (*myPositions)[aout[3]]);
          }
          //output sparse matrix
          for (int ii = 0; ii < 4; ii++){
            for (int kk = 0; kk < 4; kk++) {
              Matrix3By3 rhd = hi(ii, kk);
              outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]],
                  sqrtMass[aout[kk]], rhd, blocks[ar0].Rows, blocks[ar0].arrayPointer());

            }
          }
        }

      }
    }
  }
  
  //Dihedrals
  if (myDihedral) {
    HessDihedral hd;        //create dihedral hessian    
    unsigned int dihedrals_size = myTopo->dihedrals.size();
    for (unsigned int i = 0; i < dihedrals_size; i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->dihedrals[i].multiplicity; j++)
        if (myTopo->dihedrals[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = myTopo->dihedrals[i].atom1; aout[1] =
          myTopo->dihedrals[i].atom2;
        aout[2] = myTopo->dihedrals[i].atom3; aout[3] =
          myTopo->dihedrals[i].atom4;
        //test all in same block
        int ar0 = atom_block[aout[0]];
        if(atom_block[aout[1]] == ar0 && atom_block[aout[2]] == ar0 && atom_block[aout[3]] == ar0){
          if(simuMin){
            Real phi = dihedralAngle(aout, *myPositions);
            //
            Torsion currTorsion = myTopo->dihedrals[i];
            
            for (int j = 0; j < currTorsion.multiplicity; j++){
                currTorsion.phaseShift[j] = M_PI - currTorsion.periodicity[j] * phi;
            }
            hd.evaluate(currTorsion, (*myPositions)[aout[0]],
                        (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                        (*myPositions)[aout[3]]);
          }else{
            //
            hd.evaluate(myTopo->dihedrals[i], (*myPositions)[aout[0]],
                        (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                        (*myPositions)[aout[3]]);
          }
          //output sparse matrix
          for (int ii = 0; ii < 4; ii++){
            for (int kk = 0; kk < 4; kk++) {
              Matrix3By3 rhd = hd(ii, kk);
              outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]],
                  sqrtMass[aout[kk]], rhd, blocks[ar0].Rows, blocks[ar0].arrayPointer());
            }
          }
        }

      }
    }
  }

  //Bonds
  if (myBond){    
    unsigned int bonds_size = myTopo->bonds.size();
    for (i = 0; i < bonds_size; i++) {
      a1 = myTopo->bonds[i].atom1; a2 = myTopo->bonds[i].atom2;
      //test all in same block
      int ar1 = atom_block[a1];
      if(atom_block[a1] == atom_block[a2]){
        Real r_0 = myTopo->bonds[i].restLength;
        //
        if(simuMin) r_0 = ((*myPositions)[a2] - (*myPositions)[a1]).norm();
        //
        Real k = myTopo->bonds[i].springConstant;

        Matrix3By3 bondHess12 =
          reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
        //
        //output
        int aout[2]={a1,a2};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            if(ii == kk) rha = bondHess12;
            else rha = -bondHess12;
            outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]], 
                          sqrtMass[aout[kk]], rha, blocks[ar1].Rows, blocks[ar1].arrayPointer());
          }
        }
      }

    }
  }

  //Angles
  if (myAngle){    
    unsigned int angles_size = myTopo->angles.size();
    for (i = 0; i < angles_size; i++) {
      a1 = myTopo->angles[i].atom1;
      a2 = myTopo->angles[i].atom2;
      a3 = myTopo->angles[i].atom3;
      //test all in same block
      int ar1 = atom_block[a1];
      int aout[3];
      aout[0] = a1; aout[1] = a2; aout[2] = a3;
      //
      if(atom_block[a2] == ar1 && atom_block[a3] == ar1){
        Real theta0 = myTopo->angles[i].restAngle;
        Real k_t = myTopo->angles[i].forceConstant;
        Real ubConst = myTopo->angles[i].ureyBradleyConstant;
        Real ubRestL = myTopo->angles[i].ureyBradleyRestLength;
        //
        if(simuMin){
            Vector3D rij((*myPositions)[a2] - (*myPositions)[a1]);
            Vector3D rkj((*myPositions)[a2] - (*myPositions)[a3]);
            theta0 = atan2((rij.cross(rkj)).norm(), rij.dot(rkj));
            ubRestL = ((*myPositions)[a3] - (*myPositions)[a1]).norm();
        }
        //
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
        for (int ii = 0; ii < 3; ii++){
          for (int kk = 0; kk < 3; kk++) {
            rha = rh(ii, kk);
            outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]],
                sqrtMass[aout[kk]], rha, blocks[ar1].Rows, blocks[ar1].arrayPointer());
          }
        }
      }

    }
  }
  
#ifdef BLOCKPAIRWISEINTERACTION
  
  unsigned int _N = myTopo->atoms.size();
  //Pairwise intra block or adjacent
  
#ifdef BLOCKSCPISM
  //SCPISM
  //Pre-calculate Born radii if Self energy Hessian required
  if(myBornRadii && myBornSelf && myTopo->doSCPISM)
    evaluateBornRadii(myPositions, myTopo);
#endif
  
#ifdef BLOCKGB
  //Pre calculate Born radii for GB if required
  if (myGBBornBurialTerm && myTopo->doGBSAOpenMM) {
    evaluateGBBornBurialTerm(myPositions, myTopo);
  }
  
  if (myGBBornRadii && myTopo->doGBSAOpenMM) {
    evaluateGBBornRadii(myPositions, myTopo);
  }
#endif
  
  for (unsigned int i = 0; i < _N; i++){
    for (unsigned int j = i + 1; j < _N; j++){             
      if(atom_block[i] == atom_block[j]){  //same block}
        
        int ar1 = atom_block[i];
        
        Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        //Accumulate pairwise
        if (myLennardJones)   //Lennard jones
          rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, true);
        if (myCoulomb)        //Coulombic
          rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, true);
        
#ifdef BLOCKSCPISM
        if (myCoulombDielec)  //Coulombic Implicit solvent
          rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, true);
        if (myCoulombSCPISM)  //SCP
          rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, true);
        
        if (myBornRadii && myBornSelf && myTopo->doSCPISM)  //Bourn radii         
          rhp += evaluateBornSelfPair(i, j, myPositions, myTopo);
#endif
        
#ifdef BLOCKGB
        //GB energies
        if (myGBBornBurialTerm && myGBBornRadii && myGBACEForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBACEPair(i, j, myPositions, myTopo);
        }
        
        if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBPair(i, j, myPositions, myTopo);
        }
#endif
        
        //Output matrix
        int aout[2]={i,j};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            Matrix3By3 rhb;
            if(ii == kk) rhb = rhp;
            else rhb = -rhp;
            //outputBlocks(aout[ii], aout[kk], rhb);
            outputMatrix(atom_block_num[aout[ii]], atom_block_num[aout[kk]], sqrtMass[aout[ii]],
                         sqrtMass[aout[kk]], rhb, blocks[ar1].Rows, blocks[ar1].arrayPointer());
          }
        }
      }
    }
  }
#endif
  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate Dihedral/Improper angles
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Real BlockHessian::dihedralAngle(const int *aout, 
                                   const Vector3DBlock &myPositions) const{
  Vector3D r12(myPositions[aout[1]] - myPositions[aout[0]]);
  Vector3D r23(myPositions[aout[2]] - myPositions[aout[1]]);
  Vector3D r34(myPositions[aout[3]] - myPositions[aout[2]]);
  Vector3D a(r12.cross(r23));
  Vector3D b(r23.cross(r34));
  Vector3D c(r23.cross(a));
  a /= a.norm();
  b /= b.norm();
  c /= c.norm();
  Real cosPhi = a.dot(b);
  Real sinPhi = c.dot(b);
  return atan2(sinPhi, cosPhi);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Clear block Hessians
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::clearBlocks() {
  //Main diagonal and adjacent blocks
  for(int i=0;i<num_blocks;i++){
    blocks[i].clear();  //clear blocks
    if(i < num_blocks - 1) adj_blocks[i].clear();  //clear adjacent block
  }
  non_adj_bond_blocks.resize(0);  //non adjacent bond blocks
  non_adj_bond_index.resize(0);
  //Off diagonal bonds OR fulle electrostatic Hessian
  if(!fullElectrostatics){ //if cuttoff small do full calculation
    adj_nonbond_blocks.resize(0);  
    adj_nonbond_index.resize(0);
  }else{
    electroStatics.clear();
  }
  //memory diags
  memory_blocks = 0;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Evaluate intra and adjacent block Hessians
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::evaluateBlocks(const Real cutoffDistance, const Vector3DBlock *myPositions,
                       GenericTopology *myTopo) {
  int a1, a2, a3;
  Matrix3By3 rha;

  //
  sz = 3 * myPositions->size();
  //Impropers
  if (myImproper) outputTorsions(myTopo->impropers, *myPositions);

  //Dihedrals
  if (myDihedral) outputTorsions(myTopo->dihedrals, *myPositions);

  //RBDihedrals
  if (myRBDihedral){
    
    HessDihedral hd;        //create dihedral hessian
    //loop over all
    for (unsigned int i = 0; i < myTopo->rb_dihedrals.size(); i++) {
      RBTorsion rbt = myTopo->rb_dihedrals[i];
      bool nonZForce = false;       //test for force constants
      
      if (rbt.C0 || rbt.C1 || rbt.C2 || rbt.C3 || rbt.C4 || rbt.C5) nonZForce = true;
      
      if (nonZForce) {
        int aout[4];
        aout[0] = rbt.atom1; aout[1] = rbt.atom2;
        aout[2] = rbt.atom3; aout[3] = rbt.atom4;
        
        //test all in same block
        //
        hd.evaluate(rbt, (*myPositions)[aout[0]],
                        (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                        (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++){
          for (int kk = 0; kk < 4; kk++) {
            Matrix3By3 rhd = hd(ii, kk);
            outputBlocks(aout[ii], aout[kk], rhd);
          }
        }
                  
      }
      
    }
    
  }

  //Bonds
  if (myBond){    
    int bonds_size = myTopo->bonds.size();

    for (int i = 0; i < bonds_size; i++) {
      a1 = myTopo->bonds[i].atom1; a2 = myTopo->bonds[i].atom2;
      Real r_0 = myTopo->bonds[i].restLength;
      Real k = myTopo->bonds[i].springConstant;
      Matrix3By3 bondHess12 =
        reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
      //output sparse matrix
      //if(abs(ar0 - ar1)>=2) report << hint << " ar0 "<<ar0<<" ar1 "<<ar1<<" b1 "<<block_num_1<<" (b1) "<<atom_block_num[a1]<<" b2 "<<block_num_2<<" (b2) "<<atom_block_num[a2]<<" b1 max "<<blocks_max[ar0]<<" b2 max "<<blocks_max[ar1]<<endr;
      if(abs(atom_block[a2] - atom_block[a1])<2){
        int aout[2]={a1,a2};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            if(ii == kk) rha = bondHess12;
            else rha = -bondHess12;            
            outputBlocks(aout[ii], aout[kk], rha);
          }
        }
      }else{  
        //non-adjacent block bonds!
        int ar0 = atom_block[a1]; int ar1 = atom_block[a2];
        int block_num_1, block_num_2;
        if(ar0 > ar1){
          int tempi = ar0;
          ar0 = ar1;
          ar1 = tempi;
        }
        if(atom_block[a1] == ar0){
          block_num_1 = atom_block_num[a1];
          block_num_2 = atom_block_num[a2];
        }else{
          block_num_1 = atom_block_num[a2];
          block_num_2 = atom_block_num[a1];
        }
        BlockMatrix M((hess_eig_point[ar0]+block_num_1)*3,(hess_eig_point[ar1]+block_num_2)*3,3,3);
        M.clear();
        memory_blocks += 9;
        //output
        int aout[2]={a1,a2};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            if(ii == kk) rha = bondHess12;
            else rha = -bondHess12;
            int block_num_1 = atom_block_num[aout[ii]]; int block_num_2 = atom_block_num[aout[kk]];
            int ar0 = atom_block[aout[ii]]; int ar1 = atom_block[aout[kk]];
            if(ar0 == ar1) outputMatrix(block_num_1, block_num_2, sqrtMass[aout[ii]],
                                              sqrtMass[aout[kk]], rha, blocks[atom_block[aout[ii]]].Rows, blocks[atom_block[aout[ii]]].arrayPointer());
            if(ar0 < ar1) outputMatrix(0, 0, sqrtMass[aout[ii]],
                                              sqrtMass[aout[kk]], rha, M.Rows, M.arrayPointer());
          }
        }
        non_adj_bond_blocks.push_back(M);
        non_adj_bond_index.push_back(ar0);
        non_adj_bond_index.push_back(ar1);
      }

    }
  }

  //Angles
  if (myAngle) { 
    ReducedHessAngle rh;
    int angles_size = myTopo->angles.size();
    for (int i = 0; i < angles_size; i++) {
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
      //Find blocks
      int aout[3] = {a1, a2, a3};
      //output sparse matrix
      for (int ii = 0; ii < 3; ii++)
        for (int kk = 0; kk < 3; kk++) {
          rha = rh(ii, kk);
          outputBlocks(aout[ii], aout[kk], rha);
        }
    }
  }
  
  unsigned int _N = myTopo->atoms.size();
  //Pairwise intra block or adjacent
  
  //SCPISM
  //Pre-calculate Born radii if Self energy Hessian required
  if(myBornRadii && myBornSelf && myTopo->doSCPISM)
    evaluateBornRadii(myPositions, myTopo);
  
#ifdef ADDGB
  //Pre calculate Born radii for GB if required
  if (myGBBornBurialTerm && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Calculate forces GBBornBurialTerm"<<endr;
    evaluateGBBornBurialTerm(myPositions, myTopo);
    
  }
  
  if (myGBBornRadii && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Calculate forces GBBornRadii"<<endr;
    evaluateGBBornRadii(myPositions, myTopo);
  }
  
  //if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Appropriate flags set for calculation of GB hessian"<<endr;
  //}
  
#endif

  for (unsigned int i = 0; i < _N; i++){
    for (unsigned int j = i + 1; j < _N; j++){             
      if(abs(atom_block[i] - atom_block[j]) < 2){  //within block or adjacent
        Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        //Accumulate pairwise
        if (myLennardJones)   //Lennard jones
          rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, true);
        if (myCoulomb)        //Coulombic
          rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, true);
        if (myCoulombDielec)  //Coulombic Implicit solvent
          rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, true);
        if (myCoulombSCPISM)  //SCP
          rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, true);
#ifdef ADDSCPISM

        if (myBornRadii && myBornSelf && myTopo->doSCPISM)  //Bourn radii         
          rhp += evaluateBornSelfPair(i, j, myPositions, myTopo);
#endif
        
#ifdef ADDGB
        //GB energies
        if (myGBBornBurialTerm && myGBBornRadii && myGBACEForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBACEPair(i, j, myPositions, myTopo);
        }
        
        if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBPair(i, j, myPositions, myTopo);
        }
        
#endif

        //Output matrix
        int aout[2]={i,j};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            Matrix3By3 rhb;
            if(ii == kk) rhb = rhp;
            else rhb = -rhp;
            outputBlocks(aout[ii], aout[kk], rhb);
          }
        }
      }
    }
  }
  //
  //Do adjacent non-bonded blocks, with cutoff and by residue.  
  if(!fullElectrostatics){ //if cuttoff small do full calculation
    //adj_nonbond_blocks.resize(0);  
    //adj_nonbond_index.resize(0);
    //coarse
    int bCount = 0;
    for(int res_a=0;res_a<num_residues;res_a++){
      for(int res_b=res_a+2;res_b<num_residues;res_b++){
        int rac_a = residues_alpha_c[res_a]; int rac_b = residues_alpha_c[res_b];
        Real ac_dist;
        if(abs(atom_block[rac_a] - atom_block[rac_b]) >= 2 &&
          (ac_dist = ((*myPositions)[rac_a] - (*myPositions)[rac_b]).norm()) < cutoffDistance){
          bCount++;
          report << debug(3) << "Residues " << res_a << ", " << res_b << ", distance " << ac_dist << "." << endl;
          //
          int a0 = residues[res_a*MAX_ATOMS_PER_RES];
          int a1 = residues[res_b*MAX_ATOMS_PER_RES];
          int ar0 = atom_block[a0];
          int ar1 = atom_block[a1];
          int block_num_1 = atom_block_num[a0];
          int block_num_2 = atom_block_num[a1];
          int r_max_a3 = residues_max[res_a]*3; 
          int r_max_b3 = residues_max[res_b]*3; 
          BlockMatrix M((hess_eig_point[ar0]+block_num_1)*3,(hess_eig_point[ar1]+block_num_2)*3,r_max_a3,r_max_b3);
          M.clear();
          memory_blocks += r_max_a3 * r_max_b3;
          //
          for(int c=0;c<residues_max[res_a];c++){
            int i = residues[res_a*MAX_ATOMS_PER_RES+c];
            for(int d=0;d<residues_max[res_b];d++){
              int j = residues[res_b*MAX_ATOMS_PER_RES+d];
              //          
              Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
              if (myLennardJones)   //Lennard jones
                rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, true);
              if (myCoulomb)        //Coulombic
                rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, true);
              if (myCoulombDielec)  //Coulombic Implicit solvent
                rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, true);
              if (myCoulombSCPISM)  //SCP
                rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, true);
#ifdef ADDSCPISM
              if (myBornRadii && myBornSelf && myTopo->doSCPISM)  //Bourn radii         
                rhp += evaluateBornSelfPair(i, j, myPositions, myTopo);
#endif
              
#ifdef ADDGB
              //add GB forces
              if (myGBBornBurialTerm && myGBBornRadii && myGBACEForce && myTopo->doGBSAOpenMM) {
                rhp += evaluateGBACEPair(i, j, myPositions, myTopo);
              }
              
              if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
                rhp += evaluateGBPair(i, j, myPositions, myTopo);
              }
              
#endif

              //Output matrix
              int aout[2]={i,j};
              for (int ii = 0; ii < 2; ii++){
                for (int kk = 0; kk < 2; kk++) {
                  Matrix3By3 rhb;
                  if(ii == kk) rhb = rhp;
                  else rhb = -rhp;
                  int block_num_1 = atom_block_num[aout[ii]]; int block_num_2 = atom_block_num[aout[kk]];
                  int ar0 = atom_block[aout[ii]]; int ar1 = atom_block[aout[kk]];
                  if(ar0 == ar1) outputMatrix(block_num_1, block_num_2, sqrtMass[aout[ii]],
                                                    sqrtMass[aout[kk]], rhb, blocks[atom_block[aout[ii]]].Rows, blocks[atom_block[aout[ii]]].arrayPointer());
                  if(ar0 < ar1) outputMatrix(atom_res_num[aout[ii]], atom_res_num[aout[kk]], sqrtMass[aout[ii]],
                                                    sqrtMass[aout[kk]], rhb, M.Rows, M.arrayPointer());
                }
              }
            }
          }
          //
          adj_nonbond_blocks.push_back(M);
          adj_nonbond_index.push_back(ar0);
          adj_nonbond_index.push_back(ar1);
          //
        }
      }
    }
    report << debug(2) << " Total blocks " << bCount <<"."<<endl;
  }else{
    //Do full calculation if set
    evaluateInterBlocks(myPositions, myTopo);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Torsion outputs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::outputTorsions(const std::vector<Torsion> &torsions, const Vector3DBlock &myPositions){    
  
  HessDihedral hi;        //create improper hessian
  unsigned int torsions_size = torsions.size();
  for (unsigned int i = 0; i < torsions_size; i++) {
    bool nonZForce = false;       //test for force constants
    for (int j = 0; j < torsions[i].multiplicity; j++)
      if (torsions[i].forceConstant[j]) nonZForce = true;

    if (nonZForce) {
      int aout[4] = {torsions[i].atom1, torsions[i].atom2,
                      torsions[i].atom3, torsions[i].atom4};
      //
      hi.evaluate(torsions[i], myPositions[aout[0]],
                  myPositions[aout[1]], myPositions[aout[2]],
                  myPositions[aout[3]]);
      //output sparse matrix
      for (int ii = 0; ii < 4; ii++)
        for (int kk = 0; kk < 4; kk++) {
          Matrix3By3 rhd = hi(ii, kk);
          outputBlocks(aout[ii], aout[kk], rhd);
        }

    }
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Output 3X3 Hessian to Hessian blocks
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::outputBlocks(const unsigned int i, const unsigned int j, const Matrix3By3 rhd){
  int ar0 = atom_block[i]; int ar1 = atom_block[j];
  if(ar0 == ar1)
    outputMatrix(atom_block_num[i], atom_block_num[j], sqrtMass[i],
                  sqrtMass[j], rhd, (blocks[ar0]).Rows, blocks[ar0].arrayPointer());
  if(ar0 < ar1)
    outputMatrix(atom_block_num[i], atom_block_num[j], sqrtMass[i],
                  sqrtMass[j], rhd, adj_blocks[ar0].Rows, adj_blocks[ar0].arrayPointer());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Output inter block Hessians
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::evaluateInterBlocks(const Vector3DBlock *myPositions,
                       GenericTopology *myTopo) {
  
  unsigned int _N = myTopo->atoms.size();
  sz = 3 * myPositions->size();
  //Set Matrix
  //electroStatics.clear();
  //
  //Pairwise  
  
  //SCPISM
  //Pre-calculate Born radii if Self energy Hessian required
  if(myBornRadii && myBornSelf && myTopo->doSCPISM)
    evaluateBornRadii(myPositions, myTopo);
  
#ifdef ADDGB
  //Pre calculate Born radii for GB if required
  if (myGBBornBurialTerm && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Calculate forces GBBornBurialTerm"<<endr;
    evaluateGBBornBurialTerm(myPositions, myTopo);
    
  }
  
  if (myGBBornRadii && myTopo->doGBSAOpenMM) {
    //report << plain <<"Hessian : Calculate forces GBBornRadii"<<endr;
    evaluateGBBornRadii(myPositions, myTopo);
  }
  
#endif
  
  for (unsigned int i = 0; i < _N; i++){
    for (unsigned int j = i + 1; j < _N; j++){             
      if(abs(atom_block[i] - atom_block[j]) >= 2){  //NOT within block or adjacent
        Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        //Accumulate pairwise
        if (myLennardJones)   //Lennard jones
          rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, true);
        if (myCoulomb)        //Coulombic
          rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, true);
        if (myCoulombDielec)  //Coulombic Implicit solvent
          rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, true);
        if (myCoulombSCPISM)  //SCP
          rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, true);
#ifdef ADDSCPISM
        if (myBornRadii && myBornSelf && myTopo->doSCPISM)  //Bourn radii         
          rhp += evaluateBornSelfPair(i, j, myPositions, myTopo);
#endif
        
#ifdef ADDGB
        // add GB forces
        if (myGBBornBurialTerm && myGBBornRadii && myGBACEForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBACEPair(i, j, myPositions, myTopo);
        }
        
        if (myGBBornBurialTerm && myGBBornRadii && myGBForce && myTopo->doGBSAOpenMM) {
          rhp += evaluateGBPair(i, j, myPositions, myTopo);
        }
        
#endif

        //Output matrix
        int aout[2]={i,j};
        for (int ii = 0; ii < 2; ii++){
          for (int kk = 0; kk < 2; kk++) {
            Matrix3By3 rhb;
            if(ii == kk) rhb = rhp;
            else rhb = -rhp;
            outputMatrix(aout[ii], aout[kk], sqrtMass[aout[ii]], sqrtMass[aout[kk]], rhb, sz, electroStatics.arrayPointer());
          }
        }
      }
    }
  }

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Put matrix into Hessian array
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::outputMatrix(int i, int j, Real sqrtMassi, Real sqrtMassj, 
                                    Matrix3By3 rha, int arrSz, double *basePoint){

    for (int ll = 0; ll < 3; ll++)
        for (int mm = 0; mm < 3; mm++)
          basePoint[(i * 3 + ll) + (j * 3 + mm) * arrSz] +=
                        ((rha(ll, mm) / (sqrtMassi * sqrtMassj)));//sqrt

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// find forces in block for 'small' numerical Hessian
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::evaluateBlockForces( const unsigned int blockStart, const unsigned int blockEnd,
                            const Vector3DBlock *myPositions,
                                const GenericTopology *myTopo, Vector3DBlock *blockForces){

  //clear array
  for (unsigned int i = blockStart; i <= blockEnd; i++){
      (*blockForces)[i] = Vector3D(0.0, 0.0, 0.0);
  }

  //get boundary and dummy energies
  const VacuumBoundaryConditions &boundary =
            ((SemiGenericTopology<VacuumBoundaryConditions> &)(*myTopo)).boundaryConditions;
    
  ScalarStructure energies;
  energies.clear();

    //~~~~Bonds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myBond){
    for (unsigned i = 0; i < myTopo->bonds.size(); i++) {
      const unsigned int a1 = myTopo->bonds[i].atom1;
      const unsigned int a2 = myTopo->bonds[i].atom2;

      //in this block?
      if( (a1 >= blockStart && a1 <= blockEnd) && (a2 >= blockStart && a2 <= blockEnd) ){

        //call force calculaton
        BondSystemForce<VacuumBoundaryConditions> bsf;

        bsf.calcBond(boundary,
                        myTopo->bonds[i], myPositions, blockForces, &energies);

      }
    }
  }

  //~~~~Angles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myAngle){
    for (unsigned i = 0; i < myTopo->angles.size(); i++) {
      const unsigned int a1 = myTopo->angles[i].atom1;
      const unsigned int a2 = myTopo->angles[i].atom2;
      const unsigned int a3 = myTopo->angles[i].atom3;
      
      //in this block?
      if( (a1 >= blockStart && a1 <= blockEnd) &&
            (a2 >= blockStart && a2 <= blockEnd) &&
                (a3 >= blockStart && a3 <= blockEnd) ){
      
          AngleSystemForce<VacuumBoundaryConditions> asf;

          asf.calcAngle(boundary, myTopo->angles[i],
                        myPositions, blockForces, &energies);
          
      }
    }
  }

  //~~~~Impropers~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myImproper) {
    for (unsigned int i = 0; i < myTopo->impropers.size(); i++) {
      const unsigned int a1 = myTopo->impropers[i].atom1;
      const unsigned int a2 = myTopo->impropers[i].atom2;
      const unsigned int a3 = myTopo->impropers[i].atom3;
      const unsigned int a4 = myTopo->impropers[i].atom4;

      //in this block?
      if( (a1 >= blockStart && a1 <= blockEnd) &&
            (a2 >= blockStart && a2 <= blockEnd) &&
                (a3 >= blockStart && a3 <= blockEnd) &&
                    (a4 >= blockStart && a4 <= blockEnd) ){
          ImproperSystemForce<VacuumBoundaryConditions> isf;

          isf.calcTorsion(boundary, myTopo->impropers[i], myPositions,
                            blockForces, (energies)[ScalarStructure::IMPROPER], &energies);
      }

    }
  }

  //~~~~Dihedrals~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myDihedral) {
    for (unsigned int i = 0; i < myTopo->dihedrals.size(); i++) {
      const unsigned int a1 = myTopo->dihedrals[i].atom1;
      const unsigned int a2 = myTopo->dihedrals[i].atom2;
      const unsigned int a3 = myTopo->dihedrals[i].atom3;
      const unsigned int a4 = myTopo->dihedrals[i].atom4;

      //in this block?
		if( (a1 >= blockStart && a1 <= blockEnd) &&
                (a2 >= blockStart && a2 <= blockEnd) &&
                    (a3 >= blockStart && a3 <= blockEnd) &&
                        (a4 >= blockStart && a4 <= blockEnd) ){
			DihedralSystemForce<VacuumBoundaryConditions> dsf;

          dsf.calcTorsion(boundary, myTopo->dihedrals[i], myPositions,
                            blockForces, (energies)[ScalarStructure::DIHEDRAL], &energies);
      }

    }
  }

  //~~~~RBDihedrals~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myRBDihedral) {
    for (unsigned int i = 0; i < myTopo->rb_dihedrals.size(); i++) {
      const unsigned int a1 = myTopo->rb_dihedrals[i].atom1;
      const unsigned int a2 = myTopo->rb_dihedrals[i].atom2;
      const unsigned int a3 = myTopo->rb_dihedrals[i].atom3;
      const unsigned int a4 = myTopo->rb_dihedrals[i].atom4;

      //in this block?
		if( (a1 >= blockStart && a1 <= blockEnd) &&
                (a2 >= blockStart && a2 <= blockEnd) &&
                    (a3 >= blockStart && a3 <= blockEnd) &&
                        (a4 >= blockStart && a4 <= blockEnd) ){
			RBDihedralSystemForce<VacuumBoundaryConditions> dsf;

          dsf.calcRBTorsion(boundary, myTopo->rb_dihedrals[i], myPositions,
                            blockForces, (energies)[ScalarStructure::OTHER], &energies);
      }

    }
  }

  //~~~~Lennard Jones~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myLennardJones) {
    //loop over all pairs in this block only
    for (unsigned int i = blockStart; i <= blockEnd; i++){
        for (unsigned int j = i + 1; j <= blockEnd; j++){

          //if not bonded/dihedral
          const ExclusionClass ec = myTopo->exclusions.check(i, j);
          
          if (ec != EXCLUSION_FULL) {

            Vector3D rij = myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);

            const Real a = rij.normSquared();

            Real rawE = 0.0; Real rawForce = 0.0;

            LennardJonesForce jlf;

            jlf(rawE, rawForce, a, 1.0 / a, rij, myTopo, i, j, ec);

            // Add this force into the atom forces.
            Vector3D fij = -rij * rawForce;
            (*blockForces)[i] += fij;
            (*blockForces)[j] -= fij;

          }

        }
    }
  }

  //~~~~Coulomb~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myCoulomb) {
    //loop over all pairs in this block only
    for (unsigned int i = blockStart; i <= blockEnd; i++){
        for (unsigned int j = i + 1; j <= blockEnd; j++){

          //if not bonded/dihedral
          const ExclusionClass ec = myTopo->exclusions.check(i, j);

          if (ec != EXCLUSION_FULL) {

            Vector3D rij = myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);

            const Real a = rij.normSquared();

            Real rawE = 0.0; Real rawForce = 0.0;

            CoulombForce cf;

            cf(rawE, rawForce, a, 1.0 / a, rij, myTopo, i, j, ec);

            // Add this force into the atom forces.
            Vector3D fij = -rij * rawForce;
            (*blockForces)[i] += fij;
            (*blockForces)[j] -= fij;

          }

        }
    }

  }

  //~~~~End~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Evaluate course high frequency Hessians
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void BlockHessian::evaluateNumericalResidues(const Vector3DBlock *myPositions,
                       const GenericTopology *myTopo) {

    unsigned int block_start = 0;

    const unsigned int size = myPositions->size();
    
    const Real epsilon = 1e-12;

    const Real inv_epsilon = 1.0 / epsilon;

    //setup vectors
    Vector3DBlock blockForces, tempPositions, blockInitialForces;

    blockForces.resize(size);
    blockInitialForces.resize(size);

    //copy positions
    tempPositions = *myPositions;

    //for each block
    for(int i=0;i<num_blocks;i++){

        //get position in Hessian
        const unsigned int rowstart = blocks[i].RowStart;
        const unsigned int colstart = blocks[i].ColumnStart;

        //clear block
        blocks[i].clear();

        //find atoms in block
        const unsigned int block_max = blocks_max[i];

        //get initial forces for block
        evaluateBlockForces( block_start, block_start + block_max - 1,
                                myPositions, myTopo, &blockInitialForces);

        //calculate forces for each atom change
        for( unsigned j=0; j<block_max; j++ ){
            //x,y,z
            for( unsigned k=0; k<3; k++ ){

                //peturb
                tempPositions[block_start + j][k] += epsilon;

                //coulumn of data force calc
                evaluateBlockForces( block_start, block_start + block_max - 1,
                                        &tempPositions, myTopo, &blockForces);

                //un-peturb
                tempPositions[block_start + j][k] -= epsilon;

                //copy force array to column of numerical Hessian
                for( unsigned l=0; l<block_max; l++ ){
                    //x,y,z
                    for( unsigned m=0; m<3; m++ ){

                        blocks[i](rowstart + j*3+k, colstart + l*3+m) =
                                            -( blockForces[block_start + l][m]
                                                - blockInitialForces[block_start + l][m])
                                                    * inv_epsilon
                                                        * ( 1.0 /
                                                            ( sqrt(myTopo->atoms[block_start+l].scaledMass) * sqrt(myTopo->atoms[block_start+j].scaledMass ) )  );

                    }
                }
                
            }

        }


        //end of loop, update for next block
        block_start += block_max;

    }

}
