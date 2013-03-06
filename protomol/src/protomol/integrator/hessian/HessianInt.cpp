#include <protomol/integrator/hessian/HessianInt.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/base/PMConstants.h>
#include <iostream>
#include <stdio.h>

#ifdef HAVE_LIBFAH
    #include <cbang/os/File.h>
#else
    #include <fstream>
#endif

#include <protomol/base/Lapack.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ HessianInt

const string HessianInt::keyword("hessianint");

HessianInt::HessianInt() :
  STSIntegrator() {
  eigVec = 0;
}

HessianInt::HessianInt(Real timestep, string evec_s, string eval_s,
                       string hess_s, bool sorta, int fm, bool tef,
                       bool fdi, Real evt, int bvc, int rpb, Real bct,
                       bool masswt, bool bnm, bool aparm, bool geo, bool num, Real eps,
                       ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), evecfile(evec_s),
  evalfile(eval_s), hessfile(hess_s), sortOnAbs(sorta), numberOfModes(fm),
  textEig(tef), fullDiag(fdi),
  massWeight(masswt), noseMass(bnm), autoParmeters(aparm),
  geometricfdof(geo), numerichessians(num),
  eigenValueThresh(evt), blockCutoffDistance(bct),
  blockVectorCols(bvc), residuesPerBlock(rpb), epsilon(eps) {
  eigVec = 0;
  //
  hsn.findForces(overloadedForces);         //find forces and parameters
}

HessianInt::~HessianInt() {
  int info;

  info = 0;
  int numModes = numberOfModes;
  //fix summed Hessian to average
  if (hsn.hessM != 0 && totStep > 1 && fullDiag)
    for (unsigned int i = 0; i < sz * sz; i++)
      hsn.hessM[i] /= (double)totStep;

  //
  if (eigVec != 0 && totStep &&
      (evecfile != "" || evalfile != "" || hessfile != "")) {

    if(evecfile != "" || evalfile != ""){
      if(fullDiag){
        if(hsn.hessM != 0) {
          //Full diagonalize
          report << hint << "[HessianInt::run] diagonalizing Hessian." << endr;
          //info = diagHessian(eigVec, eigVal);
          int numFound;
          blockDiag.rediagTime.start();
          info = blockDiag.diagHessian(eigVec, blockDiag.eigVal, hsn.hessM, sz, numFound);
          if (info == 0) {
            int numneg;
            for (numneg = 0; numneg < (int)sz; numneg++)
              if (blockDiag.eigVal[numneg] > 0.0) break;

            //if(numneg) numneg--;
            report << hint << "[HessianInt::run] diagonalized! Number of negative "
              "eigenvalues = " << numneg << "." << endr;
            for (unsigned int i = 0; i < sz; i++) blockDiag.eigIndx[i] = i;
            //if (sortOnAbs) absSort();
            if (sortOnAbs) blockDiag.absSort(eigVec, blockDiag.eigVal, blockDiag.eigIndx, sz);
            blockDiag.rediagTime.stop();
            //output nose value if eigvals available
            if (noseMass) report << hint << "[HessianInt::run] Nose Mass, Q = " <<
                                    calcQ() << "." << endr;
          }else{
              report << hint << "[HessianInt::run] error info = " << info << endr;
          }
        }
      }else{
        //coarse diagonalize AND calculate course Hessian
        numModes = min(numberOfModes,(unsigned int)((residues_total_eigs * 9) / 10));
        report << hint << "[HessianInt::run] " << numModes << " approximate eigenvectors found. Maximum eigenvalue " << max_eigenvalue << "." << endr;

      }
    }
    if(info == 0){
          //output eigenvec matrix/ eigenval vector/ Hessian
          outputDiagHess(numModes);
          report.precision(3);
          report <<plain<<"NML Timing: Hessian: "<<(blockDiag.hessianTime.getTime()).getRealTime()
              <<"[s] ("<<totStep<<" times), diagonalize: "<<(blockDiag.rediagTime.getTime()).getRealTime()<<
                  "[s]."<<endl;
          if(!fullDiag) report <<plain<<"NML Memory: Hessian: "<<(hsn.memory_base + hsn.memory_blocks) * sizeof(Real) / 1000000
              <<"[Mb], diagonalize: "<<blockDiag.memory_footprint * sizeof(Real) / 1000000<<
              "[Mb], vectors: "<<sz*numberOfModes*sizeof(double)/1000000<<"[Mb]."<<endl;

    }
  }
  if (eigVec != 0) delete[] eigVec;
}

void HessianInt::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  //check options are sensible
  if (hessfile != "" && (evecfile != "" || evalfile != ""))
    THROW("[HessianInt::initialize] Cannot output Hessian after Lapack "
          "diagonalization!");

	if( !Lapack::isEnabled() ){
      THROW("Block Hessian diagonalization requires Lapack libraries.");
	}

  //creat hessian arrays and initialize sz
  int _N = app->positions.size();
  sz = 3 * _N;
  report << hint << "[HessianInt::Find Hessian] sz=" << sz << endr;
  //automatically generate parameters?
  if(autoParmeters){
    numberOfModes = 3*(int)sqrt((float)sz);
    residuesPerBlock = (int)pow((double)_N,0.6) / 15;
    blockVectorCols = 10 + (int)sqrt((float)residuesPerBlock);
    blockCutoffDistance = hsn.cutOff;
    report << hint << "[HessianInt::initialize] Auto parameters: numberOfModes " << numberOfModes <<
                    ", residuesPerBlock " << residuesPerBlock <<
                    ", blockVectorCols " << blockVectorCols <<
                    ", blockCutoffDistance " << blockCutoffDistance << "." << endr;
  }
  //number modes correct?
  if(numberOfModes == 0 || numberOfModes > sz) numberOfModes = sz;
  //
  if(fullDiag){
    hsn.initialData(sz);
    hsn.clear();
  }else{
    //assign hessian array for residues, and clear.
    bool fullE = false;
    if(blockCutoffDistance == 0.0) fullE = true;
    hsn.initialResidueData(app->topology, residuesPerBlock, fullE);
  }
  //
  int vecSize = sz * sz;
  if(!fullDiag) vecSize = sz * numberOfModes;
  try{
    eigVec = new double[vecSize];
  }catch(bad_alloc&){
      report << error << "[HessianInt::initialize] Cannot allocate memory for "
             << "Eigenvectors." << endr;
  }
  //Initialize BlockHessianDiagonalize, pass BlockHessian if Blocks (not full diag)
  if(fullDiag){
    blockDiag.initialize(sz);
  }else{
    blockDiag.initialize(&hsn, sz, (StandardIntegrator *)this);
  }
  //
  totStep = 0;
}

void HessianInt::run(int numTimesteps) {
  if (numTimesteps < 1)
    return;

  if (totStep) {
    preStepModify();
    doHalfKickdoDrift();
    calculateForces();
    if (numTimesteps > 1) {
      for (int i = 1; i < numTimesteps; i++) {
        doKickdoDrift();
        calculateForces();
        //true for mass re-weight;
        if(fullDiag){
          blockDiag.hessianTime.start();	//time Hessian
            if(numerichessians){
                numericalHessian();
            }else{
                hsn.evaluate(&app->positions, app->topology, massWeight);
            }
          blockDiag.hessianTime.stop();	//stop timer
          totStep++;
        }

      }
    }
    else {
      //true for mass re-weight;
      if(fullDiag){
        blockDiag.hessianTime.start();	//time Hessian
          if(numerichessians){
              numericalHessian();
          }else{
              hsn.evaluate(&app->positions, app->topology, massWeight);
          }
        blockDiag.hessianTime.stop();	//stop timer
        totStep++;
      }
    }
    doHalfKick();
    postStepModify();
  } else {
    calculateForces();
    //Find current Hessian
    totStep++;
    //true for mass re-weight;
    if(fullDiag){ //Full diagonalize
      blockDiag.hessianTime.start();	//time Hessian
        if(numerichessians){
            numericalHessian();
        }else{
            hsn.evaluate(&app->positions, app->topology, massWeight);
        }
      blockDiag.hessianTime.stop();	//stop timer
    }else{        //coarse diagonalize
      max_eigenvalue = blockDiag.findEigenvectors(&app->positions, app->topology,
                                                  eigVec, sz, numberOfModes,
                                                  blockCutoffDistance, eigenValueThresh,
                                                  blockVectorCols,
						  geometricfdof, numerichessians);
      residues_total_eigs = blockDiag.residues_total_eigs;
    }
    report << hint << "[HessianInt::Find Hessian] Hessian found!" << endr;
  }
  //
}

//Nose mass calculation based on Chris Sweet's Thesis
Real HessianInt::calcQ() {
    Real Q, sumF;

    if(sz > 6){
        sumF = 0;
        for(unsigned int i=6; i< sz;i++){
            if(blockDiag.eigVal[i]) sumF += sqrt(3.0 / fabs(blockDiag.eigVal[i]));
        }
        Q = 0.6/(8.0 * (sz - 6)) * sumF*sumF;
    }else{
        Q=0.0;
    }
    return Q;
}

void HessianInt::outputDiagHess(int numModes) {
  unsigned int i;
#ifdef HAVE_LIBFAH
  cb::File myFile;
#else
    ofstream myFile;
#endif

  if (hessfile != "") {
    myFile.open(hessfile.c_str(), ofstream::out);
    myFile.precision(10);
    if(fullDiag){
      //output hessian matrix to sparse form
      for (i = 0; i < sz * sz; i++)
        if (hsn.hessM[i] != 0.0)
          myFile << i / sz + 1 << " " << i % sz + 1
                 << " " << hsn.hessM[i] << endl;

    }else{
      //Output block Hessians
      for(int ii=0;ii<hsn.num_blocks;ii++){
        int b_max = hsn.blocks_max[ii]*3;
        int start_r = hsn.blocks[ii].RowStart;
        int start_c = hsn.blocks[ii].ColumnStart;
        for(int jj=0;jj<b_max*b_max;jj++){
          myFile << start_r + (jj / b_max) + 1<< " " << start_c + (jj % b_max) + 1
                  << " " << hsn.blocks[ii].MyArray[jj] << endl;
        }
      }
    }
    //close file
    myFile.close();
  }
  if (evecfile != "") {
    //output eigenvec matrix
    myFile.open(evecfile.c_str(), ofstream::out);
    myFile.precision(10);
    //
    int vecnum = numModes;//sz;
    int vecpos = sz / 3;
    //
    if (textEig) {
      myFile << vecnum * vecpos << endl;
      myFile << "! eigenvectors from Protomol/Lapack " << endl;
      for (int i = 0; i < vecnum; i++)
        for (int k = 0; k < vecpos; k++) {
          myFile << k + 1;
          for (int j = 0; j < 3; j++){
            myFile << " " << eigVec[i * sz + k * 3 + j];
          }
          myFile << endl;
        }

    } else {
      int numrec = sz / 3;
      int32 vp = vecpos;
      int32 fm = numModes;
      double ev;
      if(fullDiag) ev = blockDiag.eigVal[blockDiag.eigIndx[(numrec - 1) * 3 + 2]];
      else ev = max_eigenvalue;
      //
      //		myFile  << "! eigenvectors from Protomol/Lapack "<< endl;
      if (ISLITTLEENDIAN) swapBytes(vp);
      myFile.write((char *)&vp, sizeof(int32));
      if (ISLITTLEENDIAN) swapBytes(fm);
      myFile.write((char *)&fm, sizeof(int32));
      if (ISLITTLEENDIAN) swapBytes(ev);
      myFile.write((char *)&ev, sizeof(double));
      for (int i = 0; i < (int)numModes; i++)
        for (int k = 0; k < vecpos; k++)
          //myFile << k + 1;
          for (int j = 0; j < 3; j++) {
            double evec = eigVec[i * sz + k * 3 + j];
            if (ISLITTLEENDIAN) swapBytes(evec);
            myFile.write((char *)&evec, sizeof(double));
          }

    }
    //close file
    myFile.close();
  }
  if (evalfile != "") {
    //output eigenval vector
    myFile.open(evalfile.c_str(), ofstream::out);
    myFile.precision(10);

    myFile << numModes << endl;
    myFile << "! eigenvalues from Protomol/Lapack " << endl;
    for (int i = 0; i < numModes; i++) {
      myFile << i + 1 << " " << blockDiag.eigVal[i] << endl;
    }

    //close file
    myFile.close();
  }
}

void HessianInt::doHalfKickdoDrift() {
  if (anyPreDriftOrNextModify()) {
    doHalfKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] += (*myForces)[i] * h * 0.5 /
                            app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void HessianInt::doKickdoDrift() {
  if (anyPreDriftOrNextModify() || anyPreStepModify() ||
      anyPostStepModify()) {
    if (anyPreStepModify() || anyPostStepModify()) {
      doHalfKick();
      postStepModify();
      preStepModify();
      doHalfKick();
    } else doKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void HessianInt::numericalHessian() {
    
    //set epsilon now a parameter
    //const Real epsilon = 1e-6;
    
    //get current position forces
    calculateForces();
    
    Vector3DBlock orgForce = *myForces;
    
    report << hint << "[HessianInt::numericalHessian] Numerical Hessian calculation." << endr;
    
    //find each column
    for(unsigned int i=0; i<sz; i++){
        
        //perturb
        (app->positions)[i/3][i%3] += epsilon;
        
        calculateForces();
        
        //Vector3DBlock firstForce = *myForces;
        
        //reset positions
        (app->positions)[i/3][i%3] -= epsilon;

        //calculateForces();

        //reset positions
        //(app->positions)[i/3][i%3] += epsilon;

        Real divconst = 1.0 / epsilon;
        
        if(massWeight){
            divconst /= sqrt(app->topology->atoms[i/3].scaledMass);
        }

        //calculate finite difference
        Vector3DBlock col = (orgForce - *myForces) * divconst;

        //set column
        bool colset = hsn.setHessianColumn( col, i, app->topology, massWeight );
        
        if(!colset) std::cout << "Column " << i << " not set!" << std::endl;
    }

}

void HessianInt::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("eigvecFile",
               Value(evecfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Eigenvector filename")));
  parameters.push_back
    (Parameter("eigvalFile",
               Value(evalfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Eigenvalue filename")));
  parameters.push_back
    (Parameter("hessianFile",
               Value(hessfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Hessian sparse filename")));
  parameters.push_back
    (Parameter("sortByAbs",
               Value(sortOnAbs, ConstraintValueType::NoConstraints()), true,
               Text("Sort eigenvalues/vectors by absolute magnitude")));
  parameters.push_back
    (Parameter("numberOfModes",
               Value(numberOfModes, ConstraintValueType::NoConstraints()),
               0, Text("Number of fixed modes")));
  parameters.push_back
    (Parameter("textEigFile",
               Value(textEig, ConstraintValueType::NoConstraints()), false,
               Text("Output eigenvectors as text file")));
  parameters.push_back
    (Parameter("fullDiag",
               Value(fullDiag, ConstraintValueType::NoConstraints()), false,
               Text("Full diagonalization?")));
  parameters.push_back
    (Parameter("eigenValueThresh",
               Value(eigenValueThresh,ConstraintValueType::NotNegative()),5.0,
               Text("'Inner' eigenvalue inclusion threshold.")));
  parameters.push_back
    (Parameter("blockVectorCols",
               Value(blockVectorCols,ConstraintValueType::NotNegative()),0,
               Text("Target number of block eigenvector columns.")));
  parameters.push_back
    (Parameter("residuesPerBlock",
               Value(residuesPerBlock,ConstraintValueType::NotNegative()),1,
               Text("Residues per block.")));
  parameters.push_back
    (Parameter("blockCutoffDistance",
               Value(blockCutoffDistance,ConstraintValueType::NotNegative()),10,
               Text("Block cutoff distance for electrostatic forces.")));
  parameters.push_back
    (Parameter("massweight",
               Value(massWeight, ConstraintValueType::NoConstraints()),
               true, Text("Mass weight Hessian?")));
  parameters.push_back
    (Parameter("noseMass",
               Value(noseMass, ConstraintValueType::NoConstraints()),
               false, Text("Calculate Nose Mass?")));
  parameters.push_back
    (Parameter("autoParameters",
               Value(autoParmeters, ConstraintValueType::NoConstraints()),
               false, Text("Automatically generate diagonalization parameters.")));
  parameters.push_back
    (Parameter("geometricfdof",
               Value(geometricfdof, ConstraintValueType::NoConstraints()),
               false, Text("Calculate fixed degrees of freedom geometrically.")));
  parameters.push_back
    (Parameter("numericHessians",
               Value(numerichessians, ConstraintValueType::NoConstraints()),
               false, Text("Calculate Hessians numerically.")));
  parameters.push_back
    (Parameter("Epsilon",
               Value(epsilon, ConstraintValueType::NotNegative()),
               1e-6, Text("Epsilon for numerical Hessian.")));
}

STSIntegrator *HessianInt::doMake(const vector<Value> &values,
                                  ForceGroup *fg) const {
  return new HessianInt(values[0], values[1], values[2], values[3], values[4],
                        values[5], values[6], values[7], values[8], values[9],
                        values[10], values[11], values[12], values[13], 
                        values[14], values[15], values[16], values[17], fg);
}

