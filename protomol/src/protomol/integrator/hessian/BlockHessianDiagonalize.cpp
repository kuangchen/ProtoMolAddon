#include <protomol/integrator/hessian/BlockHessianDiagonalize.h>
#include <protomol/integrator/hessian/BlockHessian.h>

#include <protomol/base/Report.h>

#include <protomol/type/BlockMatrix.h>

#include <protomol/base/Lapack.h>

#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//Generate symetric or upper triangular Hessians
#define SYMHESS 1
//

namespace ProtoMol {
  //__________________________________________________ BlockHessianDiagonalize

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //constructors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BlockHessianDiagonalize::BlockHessianDiagonalize() {
    bHess = 0; eigVal = 0; eigIndx = 0;
    blocks_num_eigs = 0; rE = 0;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BlockHessianDiagonalize::~BlockHessianDiagonalize() {		
    if(eigVal != 0)  delete [] eigVal;
    if(eigIndx != 0) delete [] eigIndx;
    if(blocks_num_eigs != 0) delete [] blocks_num_eigs;
    if(rE != 0) delete [] rE;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize for Block Hessians
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::initialize(BlockHessian * bHessIn, 
                        const int sz, StandardIntegrator *intg) {
	this->intg = intg;
	
	if( !Lapack::isEnabled() ){
		THROW("Block Hessian diagonalization requires Lapack libraries.");
	}

    //assign pointer to BlockHessian object    
    bHess = bHessIn;  
    //assign arrays
    try{
        //assign resudue eigenvector/value array	
        blocks_num_eigs = new int[bHess->num_blocks];
        rE = new double[bHess->hess_eig_size * 3];
        //
        eigVal = new double[sz];
        eigIndx = new int[sz];
    }catch(bad_alloc&){
        report << error << "[BlockHessianDiagonalize::initialize] Block Eigenvector array allocation error." << endr;
    }
    //Assign Eigenvector Blocks
    memory_footprint = 0;
    blockEigVect.resize(bHess->num_blocks);
    for(int i=0;i<bHess->num_blocks;i++){
      unsigned int start = bHess->hess_eig_point[i]*3;
      unsigned int rows = bHess->blocks_max[i]*3;
      blockEigVect[i].initialize(start,start,rows,rows);  //initialize block
      memory_footprint += rows * rows;
    }
    //timers/counters for diagnostics
    rediagTime.reset();
    hessianTime.reset();

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize for Full Hessians
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::initialize(const int sz) {
	if( !Lapack::isEnabled() ){
		THROW("Block Hessian diagonalization requires Lapack libraries.");
	}

    //assign pointer to BlockHessian object    
    bHess = 0;  
    //assign arrays
    try{
        //
        eigVal = new double[sz];
        eigIndx = new int[sz];
    }catch(bad_alloc&){
        report << error << "[BlockHessianDiagonalize::initialize] Eigenvector array allocation error." << endr;
    }
    //timers/counters for diagnostics
    rediagTime.reset();
    hessianTime.reset();
    //
    memory_footprint = 0;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find course eigenvectors, put into array with pointer mhQu, 
  // _3N rows, _rfM columns
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Real BlockHessianDiagonalize::findEigenvectors( Vector3DBlock *myPositions,
                       GenericTopology *myTopo, 
                       double * mhQu, const int sz_row, const int sz_col, 
                       const Real blockCutoffDistance, const Real eigenValueThresh,
                       const int blockVectorCols,
                       const bool geom, const bool numeric) {

    //find 'minimum' Hessians for blocks
    hessianTime.start();	//time Hessian

    //evaluate numerically?
    if( !numeric ){
        bHess->evaluateResidues(myPositions, myTopo, false); //true for quasi-munimum, now using false as improved results
    }else{
        bHess->evaluateNumericalResidues(myPositions, myTopo);
    }

    hessianTime.stop();	//stop timer
    //
    if(OUTPUTBHESS) outputDiagnostics(2); //Output block Hessian matrix
    //Diagonalize residue Hessians       
    //find coarse block eigenvectors 
    Real max_eigenvalue = findCoarseBlockEigs(myPositions, myTopo,
                                    eigenValueThresh, blockVectorCols, geom); //true=remove geometric dof
    //
    report << debug(2) << "[BlockHessianDiagonalize::findEigenvectors] Average residue eigenvalues "<<residues_total_eigs/bHess->num_blocks<<
              ", for "<<bHess->num_blocks<<" blocks."<<endr;
    //**** Find 'inner' Hessian H *********************************************
    //Find local block Hessian, and inter block Hessian within distance 'blockCutoffDistance'
    hessianTime.start();	//time Hessian

    //~~~~inner hessian is constructed and projected here~~~~~~~~~~~~~~~~~~~~~~~

    bHess->clearBlocks();
    bHess->evaluateBlocks(blockCutoffDistance, myPositions, myTopo);
    //hessianTime.stop();	//stop timer
    
    //calculate S numerically?
    if( !numeric ){
        //put Q^T H Q into 'innerDiag' matrix
        innerHessian();
        //If cuttoff small do full calculation
        if(bHess->fullElectrostatics){
          //Q^T H Q for full electrostatics
          fullElectrostaticBlocks();
        }
    }else{
        calculateS(myPositions, myTopo);
    }

    //include inner matrix construction
    hessianTime.stop();	//stop timer

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(OUTPUTIHESS) outputDiagnostics(0); //Output inner hessian Matrices
    //**** diagonalize 'inner' Hessian ****************************************
    int numeFound;
    rediagTime.start();    
    innerEigVec.initialize(0,0,residues_total_eigs,residues_total_eigs); //set small output matrix
    int info = diagHessian(innerEigVec.arrayPointer(), 
                            eigVal, innerDiag.arrayPointer(), residues_total_eigs, numeFound);
    rediagTime.stop();
    if( info == 0 ){            
      for(int i=0;i<residues_total_eigs;i++) eigIndx[i] = i;
      absSort(innerEigVec.arrayPointer(), eigVal, eigIndx, residues_total_eigs);
    }else{
      report << error << "[BlockHessianDiagonalize::findEigenvectors] Diagnostic diagonalization failed."<<endr;
    }
    if(OUTPUTEGVAL) outputDiagnostics(1); //Output eigenvalues
    //
    //**** Multiply block eigenvectors with inner eigenvectors ****************
    innerEigVec.columnResize(sz_col);
    for(int ii=0;ii<bHess->num_blocks;ii++)
      blockEigVect[ii].productToArray(innerEigVec, mhQu, 0, 0, sz_row, sz_col); //Aaa^{T}This
    //**** Return max eigenvalue for minimizer line search ********************
    return max_eigenvalue;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find 'inner' matrix, S=Q^T H Q, using numerical methods
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::calculateS( Vector3DBlock *myPositions,
                                                GenericTopology *myTopo ){
      
    //save current time (calculate forces updates!)
    Real actTime = myTopo->time;

    //Initialize storage
    innerDiag.initialize(0,0,residues_total_eigs,residues_total_eigs); //set small output matrix
    innerDiag.clear();
    memory_footprint += residues_total_eigs * residues_total_eigs;

    //Re-size eigs space for num_eig vectors per block
    unsigned int i_res_sum = 0;
    for(int i=0;i<bHess->num_blocks;i++){
      blockEigVect[i].blockMove(blockEigVect[i].RowStart, i_res_sum);
      blockEigVect[i].columnResize(blocks_num_eigs[i]);
      i_res_sum += blocks_num_eigs[i];
    }

    //get size
    unsigned int sz = myPositions->size();

    //creat full matrix of reduced initial eigenvectors
    BlockMatrix fullEigs( 0, 0, 3 * sz, residues_total_eigs );
    fullEigs.clear();

    //copy accross
    for(int i=0;i<bHess->num_blocks;i++){
        fullEigs += blockEigVect[i];
    }
    
    //save positions
    Vector3DBlock tempPos = *myPositions;

    //get initial force
    intg->calculateForces();

    Vector3DBlock tempForce = *(intg->getForces());

    //define epsilon
    Real epsilon = 1e-9;

    //loop over each eigenvector purtubation
    for(int eg=0; eg<residues_total_eigs; eg++ ){

        //get purturbed position
        for( unsigned i=0; i < sz; i++ ){
            for( unsigned j=0; j<3; j++ ){
                (*myPositions)[i][j] = tempPos[i][j]
                    + epsilon * ( 1.0 / sqrt(myTopo->atoms[i].scaledMass))
                        * fullEigs( i*3+j, eg );

            }
        }

        //get new forces and find difference
        intg->calculateForces();

        Vector3DBlock deltaForce = *(intg->getForces()) - tempForce;

        //creat row matrix of deltaForce vector
        BlockMatrix dForce( 0, 0, 1, 3 * sz );

        for( unsigned i=0; i < 3 * sz; i++ ){
            dForce[i] = deltaForce[i/3][i%3]
                            * ( 1.0 / sqrt(myTopo->atoms[i/3].scaledMass) )
                                * ( 1.0 / epsilon );
        }

        //create output column matrix
        BlockMatrix opVector( 0, 0, 1, residues_total_eigs );
        opVector.clear();

        //get product
        for(int i=0;i<bHess->num_blocks;i++){
            dForce.product( blockEigVect[i], opVector );
        }

        //copy product
        for(int i=0;i<residues_total_eigs;i++){
            innerDiag(i,eg) = opVector[i];
        }

    }//~~~~

    //restore positions
    *myPositions = tempPos;

    //reset time
    myTopo->time = actTime;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find 'inner' matrix, S=Q^T H Q, using block matrix arithmetic
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::innerHessian(){

    //Initialize storage
    innerDiag.initialize(0,0,residues_total_eigs,residues_total_eigs); //set small output matrix
    innerDiag.clear();
    memory_footprint += residues_total_eigs * residues_total_eigs;

    //Re-size eigs space for num_eig vectors per block
    unsigned int i_res_sum = 0;
    for(int i=0;i<bHess->num_blocks;i++){
      blockEigVect[i].blockMove(blockEigVect[i].RowStart, i_res_sum);
      blockEigVect[i].columnResize(blocks_num_eigs[i]);
      i_res_sum += blocks_num_eigs[i];
    }

    //Do blocks if Hessian distance is 0
    for(int ii=0;ii<bHess->num_blocks;ii++){  
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->blocks[ii]).ColumnStart, (blockEigVect[ii]).Columns, (bHess->blocks[ii]).Columns);
      (blockEigVect[ii]).transposeProduct(bHess->blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ii], innerDiag); //Aaa^{T}HAaa
    }

    //Do blocks if Hessian distance is 1
    for(int ii=0;ii<bHess->num_blocks-1;ii++){  
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->adj_blocks[ii]).ColumnStart, (blockEigVect[ii]).Columns, (bHess->adj_blocks[ii]).Columns);
      (blockEigVect[ii]).transposeProduct(bHess->adj_blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ii+1], innerDiag); //Aaa^{T}HAbb
      if(SYMHESS){
        //Dont need this except for symmetric Hessians
        for(unsigned jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
          for(unsigned kk=(blockEigVect[ii+1]).ColumnStart;kk<(blockEigVect[ii+1]).Columns+(blockEigVect[ii+1]).ColumnStart;kk++)
             innerDiag(kk,jj) = innerDiag(jj,kk);
      }
    }
    //Do non-adjacent bond blocks, distance > 1
    int non_adj_bond_blocks_size = bHess->non_adj_bond_blocks.size();
    for(int ii=0;ii<non_adj_bond_blocks_size;ii++){  
      int ar0 = bHess->non_adj_bond_index[ii*2]; int ar1 = bHess->non_adj_bond_index[ii*2+1];
      BlockMatrix tempM((blockEigVect[ar0]).ColumnStart, (bHess->non_adj_bond_blocks[ii]).ColumnStart, (blockEigVect[ar0]).Columns, (bHess->non_adj_bond_blocks[ii]).Columns);
      (blockEigVect[ar0]).transposeProduct(bHess->non_adj_bond_blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ar1], innerDiag); //Aaa^{T}HAbb
      if(SYMHESS){
        //Dont need this except for symmetric Hessians
        for(unsigned jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
          for(unsigned kk=(blockEigVect[ar1]).ColumnStart;kk<(blockEigVect[ar1]).Columns+(blockEigVect[ar1]).ColumnStart;kk++)
            innerDiag(kk,jj) = innerDiag(jj,kk);
      }
    }
    //Do adjacent non-bond blocks
    int adj_nonbond_blocks_size = bHess->adj_nonbond_blocks.size();
    for(int ii=0;ii<adj_nonbond_blocks_size;ii++){  
      int ar0 = bHess->adj_nonbond_index[ii*2]; int ar1 = bHess->adj_nonbond_index[ii*2+1];
      BlockMatrix tempM((blockEigVect[ar0]).ColumnStart, (bHess->adj_nonbond_blocks[ii]).ColumnStart, (blockEigVect[ar0]).Columns, (bHess->adj_nonbond_blocks[ii]).Columns);
      (blockEigVect[ar0]).transposeProduct(bHess->adj_nonbond_blocks[ii], tempM); //Aaa^{T}This
      tempM.sumProduct(blockEigVect[ar1], innerDiag); //Aaa^{T}HAbb
      if(SYMHESS){
        //Dont need this except for symmetric Hessians
        for(unsigned jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
          for(unsigned kk=(blockEigVect[ar1]).ColumnStart;kk<(blockEigVect[ar1]).Columns+(blockEigVect[ar1]).ColumnStart;kk++)
            innerDiag(kk,jj) = innerDiag(jj,kk);
      }
    }
    //
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find 'inner' matrix, S'=Q^T H' Q, for full electrostatic Hessian H'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::fullElectrostaticBlocks(){
  
    //Full electrostatics
    for(int ii=0;ii<bHess->num_blocks;ii++){
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->electroStatics).ColumnStart, (blockEigVect[ii]).Columns, bHess->electroStatics.Columns);
      (blockEigVect[ii]).transposeProduct(bHess->electroStatics, tempM); //Aaa^{T}This
      int llStart = ii;
      if(SYMHESS) llStart = 0;
      for(int ll=llStart;ll<bHess->num_blocks;ll++){// should be from ii unless symmetric Hessians
        tempM.sumProduct(blockEigVect[ll], innerDiag); //Aaa^{T}HAbb
      }
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find isolated block eigenvectors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Real BlockHessianDiagonalize::findCoarseBlockEigs(const Vector3DBlock *myPositions,
                                                        const GenericTopology *myTopo,
                                                            const Real eigenValueThresh,
                                                                const int blockVectorCols,
                                                                    const bool geom) {

    //clear number of cols vector
    blocVectCol.clear();
    //number of cols precidence?
    bool blockVprec = false;
    if(blockVectorCols > 0) blockVprec = true;

    //size and position eigvects
    for(int i=0;i<bHess->num_blocks;i++){
      unsigned int start = blockEigVect[i].RowStart;
      blockEigVect[i].blockMove(start, start);  //make square
      blockEigVect[i].columnResize(blockEigVect[i].Rows);
    }
    //    
    residues_total_eigs = 0;	//find total residue modes to use
    Real max_eigenvalue = 0;
    int bHess_num_blks = bHess->num_blocks;

    //find atom block start
    unsigned int block_start = 0;

    //for each block
    for(int ii=0;ii<bHess_num_blks;ii++){
      int numFound;

      //diagonalize block
      rediagTime.start();
      int infor = diagHessian(blockEigVect[ii].arrayPointer(), &rE[bHess->hess_eig_point[ii] * 3], 
                                bHess->blocks[ii].arrayPointer(), bHess->blocks[ii].Rows, numFound);
      rediagTime.stop();

      if(infor) report << error << "[BlockHessianDiagonalize::findCoarseBlockEigs] Residue "<<ii+1<<" diagonalization failed."<<endr;

      //sort by magnitude of eigenvalue
      for(int i=0;i<bHess->blocks_max[ii] * 3;i++) eigIndx[i] = i;
      absSort(blockEigVect[ii].arrayPointer(), &rE[bHess->hess_eig_point[ii] * 3], eigIndx, bHess->blocks_max[ii] * 3);

      //~~~~Use geometrically generated conserved dof to generate a new basis set?
      if( geom ){
          //set conserved dof
          unsigned int cdof = 6;

          //create temporary matrix
          const unsigned int rowstart = blockEigVect[ii].RowStart;
          const unsigned int colstart = blockEigVect[ii].ColumnStart;

          BlockMatrix tmpEigs( rowstart, colstart,
                                blockEigVect[ii].Rows, blockEigVect[ii].Columns );

          //clear it
          tmpEigs.clear();

          //find positions
          const unsigned int block_max = bHess->blocks_max[ii];

          //find center and norm (1/sqrt mass)
          Vector3D pos_center(0.0,0.0,0.0);
          Real totalmass = 0.0;
          Real inorm = 0;

          for( unsigned jj=0; jj<block_max; jj++ ){
              const unsigned int atomindex = block_start + jj;
              //sums
              const Real mass = myTopo->atoms[atomindex].scaledMass;
              pos_center += (*myPositions)[atomindex] * mass;
              totalmass += mass;

              const Real oneosrm = sqrt(myTopo->atoms[atomindex].scaledMass);
              inorm += oneosrm * oneosrm;
          }

          //get inverse of actual norm
          inorm = 1.0 / sqrt(inorm);

          //actual center
          pos_center /= totalmass;//(Real)block_max;

          //create fixed dof vectors
          for( unsigned jj=0; jj<block_max; jj++ ){
              const unsigned int atomindex = block_start + jj;

              //translational dof
              Real factor = inorm * sqrt(myTopo->atoms[atomindex].scaledMass);
              tmpEigs(rowstart + jj*3,colstart) = factor;
              tmpEigs(rowstart + jj*3+1,colstart + 1) = factor;
              tmpEigs(rowstart + jj*3+2,colstart + 2) = factor;

              //rotational dof
              //cross product of rotation axis and vector to center of molecule
              //axb=ia2b3+ja3b1+ka1b2-ia3b2-ja1b3-ka2b1
              //x-axis (b1=1) ja3-ka2
              //y-axis (b2=1) ka1-ia3
              //z-axis (b3=1) ia2-ja1
              Vector3D diff = (*myPositions)[atomindex] - pos_center;
              //x
              tmpEigs(rowstart + jj*3+1,colstart + 3) = diff.c[2] * factor;//z;
              tmpEigs(rowstart + jj*3+2,colstart + 3) = -diff.c[1] * factor;//y;
              //y
              tmpEigs(rowstart + jj*3,colstart + 4) = -diff.c[2] * factor;//z;
              tmpEigs(rowstart + jj*3+2,colstart + 4) = diff.c[0] * factor;//x;
              //y
              tmpEigs(rowstart + jj*3,colstart + 5) = diff.c[1] * factor;//y;
              tmpEigs(rowstart + jj*3+1,colstart + 5) = -diff.c[0] * factor;//x;

          }

          //~~~~Norm vector 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          //norm the rotational vectors
          Real rotnorm = 0.0;
          for( unsigned ll=0; ll<block_max * 3; ll++ ){
                rotnorm += tmpEigs(rowstart + ll, colstart + 3)
                                * tmpEigs(rowstart + ll, colstart + 3);
          }
          rotnorm = 1.0 / sqrt(rotnorm);
          //do norm
          for( unsigned ll=0; ll<block_max * 3; ll++ ){
                tmpEigs(rowstart + ll, colstart + 3) *= rotnorm;
          }
          //~~~~orthoganalize vectors 2 and 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          for( unsigned jj=4; jj<=5; jj++ ){

              //(new vectors already have cdof vectors in them)
              for( unsigned kk=3; kk<jj; kk++ ){

                  //get dot product (and norm of dots)
                  Real dots = 0.0;
                  for( unsigned ll=0; ll<block_max * 3; ll++ ){
                    dots += tmpEigs(rowstart + ll, colstart + kk)
                                    * tmpEigs(rowstart + ll, colstart + jj);
                  }

                  //subtract it from current vector
                  for( unsigned ll=0; ll<block_max * 3; ll++ ){
                      tmpEigs(rowstart + ll,colstart + jj) -=
                              tmpEigs(rowstart + ll,colstart + kk) * dots;
                  }

              }

              Real nnorm = 0.0;
              for( unsigned ll=0; ll<block_max * 3; ll++ ){
                  nnorm += tmpEigs(rowstart + ll,colstart + jj) * tmpEigs(rowstart + ll,colstart + jj);
              }

              nnorm = 1.0 / sqrt(nnorm);
              //scale
              for( unsigned ll=0; ll<block_max * 3; ll++ ){
                  tmpEigs(rowstart + ll,colstart + jj) *= nnorm;
              }

          }
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          //sift actual vectors into array, orthogonalizing w.r.t. the conserved dof.
          //conserved dof, cdof, set at top of routine

          //loop over each vector in block
          for( unsigned jj=0; jj< 3 * block_max - cdof; jj++ ){

              //copy original vector
              for( unsigned ll=0; ll<block_max * 3; ll++ ){
                  tmpEigs(rowstart + ll, colstart + jj + cdof) = blockEigVect[ii](rowstart + ll,colstart + jj);
              }

              //get dot product with eack previous vector
              //(new vectors already have cdof vectors in them)
              for( unsigned kk=0; kk<jj+cdof; kk++ ){

                  //get dot product (and norm of dots)
                  Real dots = 0.0;
                  for( unsigned ll=0; ll<block_max * 3; ll++ ){
                    dots += tmpEigs(rowstart + ll, colstart + kk)
                                    * blockEigVect[ii](rowstart + ll, colstart + jj);
                  }

                  //subtract it from current vector
                  for( unsigned ll=0; ll<block_max * 3; ll++ ){
                      tmpEigs(rowstart + ll,colstart + jj + cdof) -=
                              tmpEigs(rowstart + ll,colstart + kk) * dots;
                  }

              }

              //normalize the residual vector
              Real nnorm = 0.0;
              for( unsigned ll=0; ll<block_max * 3; ll++ ){
                  nnorm += tmpEigs(rowstart + ll,colstart + jj + cdof) * tmpEigs(rowstart + ll,colstart + jj + cdof);
              }

              //remove if 1/20 th of original
              if( nnorm < 0.05 && cdof > 0 ){//TODO## AND not last?
                  report << debug(12) << "Residual vector " << jj << " norm is low = " << nnorm << " in block " 
                                << ii << ", skipping." << endr;
                  cdof--;
              }else{
                  
                  nnorm = 1.0 / sqrt(nnorm);
                  //scale
                  for( unsigned ll=0; ll<block_max * 3; ll++ ){
                      tmpEigs(rowstart + ll,colstart + jj + cdof) *= nnorm;
                  }
                  
              }

          }

          //copy across
          blockEigVect[ii] = tmpEigs;

          //update start atom for next block
          block_start += block_max;

      }
      //~~~~End geometric~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      //find number of eigs required, by the eigenvalue threshold, if block vector number not set
      //fix for no block vectors and eig_thresh greater than max_eig
      const unsigned int blockiisize = bHess->blocks_max[ii] * 3;
      blocks_num_eigs[ii] = blockiisize;
      
      for(unsigned jj=0;jj<blockiisize;jj++){
        
        const Real currEig = fabs(rE[bHess->hess_eig_point[ii] * 3 + jj]);
        
        if( !blockVprec ){
          if( currEig >= eigenValueThresh ){
            blocks_num_eigs[ii] = jj;
            break;
          }
        }else{
          blocVectCol.push_back(currEig);
        }
      }
      
      //fix residues_total_eigs
      residues_total_eigs += blocks_num_eigs[ii];

      //find maximum eigenvalue
      Real tempf = rE[(bHess->hess_eig_point[ii] + bHess->blocks_max[ii] - 1) * 3];
      if(max_eigenvalue < tempf) max_eigenvalue = tempf;
      //
    }

    //use target number of block eigenvectors?
    if(blockVprec && 
       (unsigned)(bHess_num_blks * blockVectorCols) < blocVectCol.size()){ 
      sort( blocVectCol.begin(), blocVectCol.end() ); //sort
      Real newEigThr = blocVectCol[bHess_num_blks * blockVectorCols];
      residues_total_eigs = 0;
      
      for(int ii=0;ii<bHess_num_blks;ii++){
        
        //find number of eigs for new thresh
        for(int jj=0;jj<bHess->blocks_max[ii] * 3;jj++){
          Real currEig = fabs(rE[bHess->hess_eig_point[ii] * 3 + jj]);
          if(currEig >= newEigThr){
            blocks_num_eigs[ii] = jj;
            residues_total_eigs += jj;
            break;
          }
        }
        
      }
      
    }


    return max_eigenvalue;
    //
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Diagnostic output
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::outputDiagnostics(int typ){    
    ofstream myFile;	//for diagnostic output

    if(typ == 0){
      //Output 'inner' matrix for diagnostics
      myFile.open("eigM", ofstream::out);
      myFile.precision(10);
      for(int jj=0;jj<residues_total_eigs*residues_total_eigs;jj++)
        myFile << jj / residues_total_eigs + 1 << " " << jj % residues_total_eigs + 1
                  << " " << innerDiag.MyArray[jj] << endl;
      myFile.close();
    }    
    if(typ == 1){
      //output eigenvalues for diagnostics
      myFile.open("eigRed", ofstream::out);
      myFile.precision(10);
      for(int jj=0;jj<residues_total_eigs;jj++)
        myFile << jj << "  " << eigVal[jj] << endl;
      //close file
      myFile.close();
    }
    if(typ == 2){
      //Output block Hessians
      myFile.open("blockH", ofstream::out);
      myFile.precision(10);    
      for(int ii=0;ii<bHess->num_blocks;ii++){
        int b_max = bHess->blocks_max[ii]*3; 
        int start_r = bHess->blocks[ii].RowStart;
        int start_c = bHess->blocks[ii].ColumnStart;
        for(int jj=0;jj<b_max*b_max;jj++){
          myFile << start_r + (jj / b_max) + 1<< " " << start_c + (jj % b_max) + 1
                  << " " << bHess->blocks[ii].MyArray[jj] << endl;
        }
      }
      myFile.close();
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Diagonalize Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int BlockHessianDiagonalize::diagHessian(double *eigVecO, double *eigValO, double *hsnhessM, int dim, int &numFound){
   double *wrkSp = 0;
   int *isuppz, *iwork;


   isuppz = new int[2*dim];
   iwork = new int[10*dim];
   //Diagonalize
//QR?
#if defined( HAVE_QRDIAG )    
   //copy Hessian to eig I/O space, could just copy the upper triangular part with diagonal
   const int dimsq = dim*dim;
   for(int i=0;i<dimsq;i++) eigVecO[i] = hsnhessM[i];
   wrkSp = new double[1];
#else
   wrkSp = new double[26*dim];
#endif
//
    char jobz = 'V'; char range = 'A'; char uplo = 'U'; /* LAPACK checks only first character N/V */
    int n = dim;             /* order of coefficient matrix a  */
    int lda = dim;           /* leading dimension of a array*/
    double vl = 1.0;
    double vu = 1.0; 
    int il = 1;
    int iu = 1;
    double abstol = 0;
    int ldz = dim; int lwork = 26*dim; /* dimension of work array*///int m; 
    int liwork = 10*dim;						/* dimension of int work array*/
    //Recomended abstol for max precision
    char cmach = 's';//String should be safe min but is shortened to remove warning
    int info = 0;				/* output 0=success */
    int m = 0;
    //call LAPACK 
    //	
//QR?
#if defined( HAVE_QRDIAG )    
    lwork = -1;
    Lapack::dsyev(&jobz, &uplo, &n, eigVecO, &lda, eigValO, wrkSp, &lwork, &info);
    if(info == 0){
      lwork = wrkSp[0];
      delete [] wrkSp;
      wrkSp = new double[lwork];
      Lapack::dsyev(&jobz, &uplo, &n, eigVecO, &lda, eigValO, wrkSp, &lwork, &info);
    }
#else
    abstol = Lapack::dlamch( &cmach);	//find machine safe minimum  
    //abstol = 1e-15;	//use small value for tolerence  
    //
    Lapack::dsyevr( &jobz, &range, &uplo, &n, hsnhessM, &lda, &vl, &vu, &il, &iu, &abstol, &m, eigValO, eigVecO, &ldz, isuppz, 
                wrkSp, &lwork, iwork, &liwork, &info);
#endif
	numFound = m;
    //delete arrays
    delete [] iwork;
    delete [] isuppz;
    delete [] wrkSp;
    //return status
    return info;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sort vectors for absolute value
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::absSort(double *eigVec, double *eigVal, int *eigIndx, int dim){
    int i;

    //find minimum abs value
    double minEv = fabs(eigVal[0]);
    for(i=1;i<dim;i++){
        if(minEv < fabs(eigVal[i])) break;
        else minEv = fabs(eigVal[i]);
    }
    i--;
    //sort around min
    if(i>0){
        int j = 0;
        eigIndx[j++] = i;
        int negp = i-1;
        int posp = i+1;
        while(negp >= 0 && posp < dim){
            if(fabs(eigVal[negp]) < fabs(eigVal[posp]))
                eigIndx[j++] = negp--;
            else eigIndx[j++] = posp++;
        }
        while(negp >= 0) eigIndx[j++] = negp--;
    }
    //Sort actual eigenvector array
    double *tmpVect = new double[dim];
    double tmpElt, tmpEval;
    int ii, k;
    for(i=0;i<dim;i++){
        if( eigIndx[i] != (int)i && eigIndx[i] != -1){		//need to swap?
            for(int j=0;j<dim;j++){
                tmpVect[j] = eigVec[i*dim+j];
                eigVec[i*dim+j] = eigVec[eigIndx[i]*dim+j];
            }
			//
			tmpEval = eigVal[i];
			eigVal[i] = eigVal[eigIndx[i]];
			//
            eigIndx[i] = -1;								//flag swapped
            ii = i;
            do{
                for(k=0;k<dim && eigIndx[k]!=(int)ii;k++);	//find where tmpVect goes
                if(k==dim || k==ii) break;					//end chain where indeces are equal
                for(int j=0;j<dim;j++){				//put it there
                    tmpElt = tmpVect[j];
                    tmpVect[j] = eigVec[k*dim+j];
                    eigVec[k*dim+j] = tmpElt;
                }
				//
				tmpElt = tmpEval;
				tmpEval = eigVal[k];
				eigVal[k] = tmpElt;
				//
                eigIndx[k] = -1;							//flag swapped
                ii = k;
            }while(k<dim);
        }
    }
    delete [] tmpVect;
  }

}
