#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/ScalarStructure.h>

#include <protomol/base/Lapack.h>

using namespace ProtoMol::Report;

namespace ProtoMol {

  //constructors
  NormalModeUtilities::NormalModeUtilities(): firstMode(1), numMode(-1), myGamma(-1), mySeed(-1), myTemp(-1)
    {/*tmpFX=NULL;*/ tmpC=NULL; invSqrtMass=NULL; sqrtMass=NULL;}

    NormalModeUtilities::NormalModeUtilities( int firstmode, int nummode, Real gamma, int seed, Real temperature):
        firstMode(firstmode), numMode(nummode), myGamma(gamma/ (1000 * Constant::INV_TIMEFACTOR)), mySeed(seed), myTemp(temperature)
  {
    /*tmpFX=NULL;*/ tmpC=NULL; invSqrtMass=NULL; sqrtMass=NULL;
  }

  NormalModeUtilities::~NormalModeUtilities()
  {
      //
      //if(tmpFX!=NULL) delete [] tmpFX;
      if(tmpC!=NULL) delete [] tmpC;
      if(invSqrtMass!=NULL) delete [] invSqrtMass;
      if(sqrtMass!=NULL) delete [] sqrtMass;
  }

  void NormalModeUtilities::initialize(int sz, ProtoMolApp *app, Vector3DBlock * myForces, int nm_flags){
	if( !Lapack::isEnabled() ){
		THROW("Normal Mode integrators require Lapack libraries.");
	}

    //Set up eigenvector pointers
    numEigvectsu = app->eigenInfo.myNumEigenvectors;
    eigVecChangedP = &app->eigenInfo.myEigVecChanged;
    eigValP = &app->eigenInfo.myMaxEigenvalue;
    Q = &app->eigenInfo.myEigenvectors;

    //find topology pointer
    GenericTopology *myTopo = app->topology;

    //local force pointer
    myForcesP = myForces;
    //type
    if(nm_flags & COMPLIMENT_FORCES) complimentForces = true;
    else  complimentForces = false;
    if(nm_flags & GEN_COMP_NOISE) genCompNoise = true;
    else  genCompNoise = false;
    //
    _N = sz;
    _3N = 3*_N;	//used everywhere, must be initialized
    //setup for auto re-diag by allowing full set of eigenvectors
    if(!numEigvectsu) numEigvectsu = _3N;
    //first mode?
    if(firstMode < 1 || firstMode > _3N) report << error << "firstmode = "<<firstMode<<", must be > 0 and less than 3N = "<<_3N<<"."<<endr;
    //check numMode
    if(numMode < 1 || numMode > _3N - (firstMode-1)) report << error << "numbermodes = "<<numMode<<", must be > 0 and less than or equal to 3N-(firstmode-1) = "<<_3N-(firstMode-1)<<"."<<endr;
    //
    if(numMode > (int)numEigvectsu - (firstMode-1))  report << error << "Insufficient eigenvectors: numbermodes = "<<numMode<<", must be less than than or equal to m-(firstmode-1) = "<<numEigvectsu-(firstMode-1)<<"."<<endr;
    //number of low frequency modes
    _rfM = numMode + (firstMode-1);
    //Degrees of freedom
    myTopo->degreesOfFreedom = _rfM - (firstMode-1);
    //Removed 11/03/2010 by CRS, correct temperature calculated when assuming NO conserved DOF: if(firstMode == 1) myTopo->degreesOfFreedom -= 3;
    //temporary mode space variable for intermediate calculations
    tmpC = new double[_3N];
    //define temporary position/force vector
    //tmpFX = new double[_3N];
    //setup sqrt masses and their inverse
    invSqrtMass = new double[_N];
    sqrtMass = new double[_N];
    for(int i=0;i<_N;i++){
        sqrtMass[i] = sqrt( myTopo->atoms[i].scaledMass );
        if(sqrtMass[i]) invSqrtMass[i] = 1.0 / sqrtMass[i];
    }
    //
    posTemp.resize(_N);
    //***********************************************************
    //gauss randoms
    gaussRandCoord1.resize(_N);
    gaussRandCoord2.resize(_N);
    if(genCompNoise) tempV3DBlk.resize(_N);	//use comliment noise gen if flagged

  }

  //test for remaining modes, equivalent to old fixModes=0
  bool NormalModeUtilities::testRemainingModes(){
    if((numMode + (firstMode-1)) < _3N) return true;
    else return false;
  }

  //Setup string of integrators so that all point to the top integrator 'Eigenpointers' :-)
  void NormalModeUtilities::setIntegratorSetPointers(Integrator *integrator, EigenvectorInfo *eipt, bool eiValid){

    NormalModeUtilities *nmint;
    if(eiValid){
        numEigvectsu = eipt->myNumEigenvectors;
    }else{
        numEigvectsu = 0;
    }
    eigVecChangedP = &eipt->myEigVecChanged;
    eigValP = &eipt->myMaxEigenvalue;
    Q = &eipt->myEigenvectors;
    //find next integrators
    for(Integrator* i = integrator->next();i != NULL;i = i->next()){
        nmint = dynamic_cast<NormalModeUtilities*>(i);
        if(nmint == NULL) report << error << "Normal Mode integrator chain contains unknown integrator type."<<endr;
        nmint->eigVecChangedP = &eipt->myEigVecChanged;
        nmint->eigValP = &eipt->myMaxEigenvalue;
        nmint->Q = &eipt->myEigenvectors;
        nmint->numEigvectsu = numEigvectsu;
    }
  }

  //*************************************************************************************
  //****Projectors for complement sub space and sub space********************************
  //*************************************************************************************

  //Find forces acting outside subspace
  Vector3DBlock* NormalModeUtilities::nonSubspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce){
    //
    if (force != iPforce) // Different vectors
      iPforce->intoAssign(*force);
    //vector3DBlockTOvect(force, tmpFX);	//put positions-x_0 into linear array
    //calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
    //f'=M^{-1/2}*f
    for( int i=0; i < _3N; i++)
         iPforce->c[i] *= invSqrtMass[i/3];
    //c=hQ^T*M^{-1/2}*f
    char transA = 'T';							// Transpose Q, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM; int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;		//multiplyers, see Blas docs.

    Lapack::dgemv(&transA, &m, &n, &alpha, (*Q), &m, iPforce->c, &incxy, &beta, tmpC, &incxy);

    //calculate f''=M^{-1/2}*f'-hQc using BLAS
    char transB = 'N';
    alpha = -1.0;	beta = 1.0;

    Lapack::dgemv(&transB, &m, &n, &alpha, (*Q), &m, tmpC, &incxy, &beta, iPforce->c, &incxy);

    //f'''=M^{1/2}*f''
    for( int i=0; i < _3N; i++)
            iPforce->c[i] *= sqrtMass[i/3];
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPforce);
    //delete temporary array
    return iPforce;
  }

  //Find positions outside subspace
  Vector3DBlock* NormalModeUtilities::nonSubspacePosition(Vector3DBlock * force, Vector3DBlock * iPforce){
    //
    if (force != iPforce)
      iPforce->intoAssign(*force);
    //vector3DBlockTOvect(force, tmpFX);	//put positions-x_0 into linear array
    //calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
    //f'=M^{-1/2}*f
    for( int i=0; i < _3N; i++)
            iPforce->c[i] *= sqrtMass[i/3];
    //c=hQ^T*M^{-1/2}*f
    char transA = 'T';							// Transpose Q, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM; int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;		//multiplyers, see Blas docs.

    Lapack::dgemv(&transA, &m, &n, &alpha, (*Q), &m, iPforce->c, &incxy, &beta, tmpC, &incxy);

    //calculate f''=M^{-1/2}*f'-hQc using BLAS
    char transB = 'N';
    alpha = -1.0;	beta = 1.0;

    Lapack::dgemv(&transB, &m, &n, &alpha, (*Q), &m, tmpC, &incxy, &beta, iPforce->c, &incxy);

    //f'''=M^{1/2}*f''
    for( int i=0; i < _3N; i++)
            iPforce->c[i] *= invSqrtMass[i/3];
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPforce);
    //delete temporary array
    return iPforce;
  }

  //Find forces acting inside subspace
  Vector3DBlock* NormalModeUtilities::subspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce){
    //
    if (force != iPforce)
      iPforce->intoAssign(*force);
    //vector3DBlockTOvect(force, tmpFX);	//put positions-x_0 into linear array
    //calculate M^{1/2}QQ^TM^{-1/2}f using BLAS
    //f'=M^{-1/2}*f
    for( int i=0; i < _3N; i++) {
            iPforce->c[i] *= invSqrtMass[i/3];
    }
    //c=Q^T*M^{-1/2}*f
    char transA = 'T';							// Transpose, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM-(firstMode-1); int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;

    Lapack::dgemv(&transA, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, iPforce->c, &incxy, &beta, tmpC, &incxy);

    //f''=Qc
    char transB = 'N'; /* LAPACK checks only first character N/V */
    alpha = 1.0;	beta = 0.0;

    Lapack::dgemv(&transB, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, tmpC, &incxy, &beta, iPforce->c, &incxy);

    //f'''=M^{1/2}*f''
    for( int i=0; i < _3N; i++) {
            iPforce->c[i] *= sqrtMass[i/3];
    }
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPforce);
    //delete temporary array
    return iPforce;
  }

  //Find velocities acting inside subspace
  Vector3DBlock* NormalModeUtilities::subspaceVelocity(Vector3DBlock * force, Vector3DBlock * iPforce){
    //
    if (force != iPforce)
      iPforce->intoAssign(*force);
    //vector3DBlockTOvect(force, tmpFX);	//put positions-x_0 into linear array
    //calculate M^{-1/2}QQ^TM^{1/2}f using BLAS
    //v'=M^{1/2}*v
    for( int i=0; i < _3N; i++)
            iPforce->c[i] *= sqrtMass[i/3];
    //c=Q^T*M^{-1/2}*v
    char transA = 'T';							//Transpose, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM-(firstMode-1); int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;

    Lapack::dgemv(&transA, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, iPforce->c, &incxy, &beta, tmpC, &incxy);

    //v''=Qc
    char transB = 'N'; /* LAPACK checks only first character N/V */
    alpha = 1.0;	beta = 0.0;

    Lapack::dgemv(&transB, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, tmpC, &incxy, &beta, iPforce->c, &incxy);

    //v'''=M^{-1/2}*v''
    for( int i=0; i < _3N; i++)
            iPforce->c[i] *= invSqrtMass[i/3];
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPforce);

    //delete temporary array
    return iPforce;
  }

  void NormalModeUtilities::subSpaceSift(Vector3DBlock *velocities, Vector3DBlock *forces){
    //sift current data into subspace
    subspaceVelocity(velocities, velocities);
    subspaceForce(forces, forces);
  }

  //force projection for post force modifier
  void NormalModeUtilities::forceProjection(){
      if((*Q) != NULL){
          if(complimentForces) nonSubspaceForce(myForcesP, myForcesP);
          else subspaceForce(myForcesP, myForcesP);
      }
  }

  //Project from mode subspace to 3D space
  Vector3DBlock* NormalModeUtilities::cartSpaceProj(double *tmpC, Vector3DBlock * iPos, Vector3DBlock * ex0){
    char transA = 'N';							// Transpose Q, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM; int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;		//multiplyers, see Blas docs.

    Lapack::dgemv(&transA, &m, &n, &alpha, (*Q), &m, tmpC, &incxy, &beta, iPos->c, &incxy);

    for( int i=0; i < _3N; i++)
            iPos->c[i] /= sqrtMass[i/3];
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPos);
    //add ex0
    for( int i=0; i < _N; i++)
        (*iPos)[i] += (*ex0)[i];
    //delete temporary array
    return iPos;
  }

  //Project from 3D space to mode subspace
  double* NormalModeUtilities::modeSpaceProj(double *cPos, Vector3DBlock * iPos, Vector3DBlock * ex0){
    //
    Vector3DBlock vPos = *iPos;
    for( int i=0; i < _N; i++) vPos[i] -= (*ex0)[i]; //subtract ex0
    //vector3DBlockTOvect(&vPos, tmpFX);	//put second positions-x_0 into linear array
    //calculate M^{-1/2}QQ^TM^{1/2}f using BLAS
    //v'=M^{1/2}*v
    for( int i=0; i < _3N; i++)
        vPos.c[i] *= sqrtMass[i/3];
    //c=Q^T*M^{-1/2}*v
    char transA = 'T';							//Transpose, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM; int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;

	Lapack::dgemv(&transA, &m, &n, &alpha, (*Q), &m, vPos.c, &incxy, &beta, cPos, &incxy);

    return cPos;
  }

  //*************************************************************************************
  //****Langevin Thermostat**************************************************************
  //*************************************************************************************

  // Generate projected vector of gausians
  void NormalModeUtilities::genProjGauss(Vector3DBlock *gaussRandCoord, GenericTopology *myTopo) {
    //generate set of random force variables and project into sub space
    for( int i = 0; i < _3N; i++ )
        (*gaussRandCoord)[i/3][i%3] = randomGaussianNumber(mySeed);//
    for( int i = 0; i < _N; i++ )
        (*gaussRandCoord)[i] *= sqrtMass[i];
    if(complimentForces) nonSubspaceForce(gaussRandCoord, gaussRandCoord);
    else subspaceForce(gaussRandCoord, gaussRandCoord);
    for( int i = 0; i < _N; i++ )
        (*gaussRandCoord)[i] /= myTopo->atoms[i].scaledMass;
  }

  // Generate projected vector of gausians AND its compliment
  void NormalModeUtilities::genProjGaussC(Vector3DBlock *gaussRandCoord, Vector3DBlock *gaussRandCoordm, GenericTopology *myTopo) {
    //generate set of random force variables and project into sub space
    if((int)gaussRandCoord->size() != _N || (int)gaussRandCoordm->size() != _N) return;
    for( int i = 0; i < _3N; i++ )
        (*gaussRandCoord)[i/3][i%3] = randomGaussianNumber(mySeed);//
    for( int i = 0; i < _N; i++ )
        (*gaussRandCoord)[i] *= sqrtMass[i];
    //get randoms for compliment
    *gaussRandCoordm = *gaussRandCoord;
    //
    if(complimentForces) nonSubspaceForce(gaussRandCoord, gaussRandCoord);
    else subspaceForce(gaussRandCoord, gaussRandCoord);
    //
    for( int i = 0; i < _N; i++ ){
        (*gaussRandCoordm)[i] -= (*gaussRandCoord)[i];
        (*gaussRandCoord)[i] /= myTopo->atoms[i].scaledMass;
        (*gaussRandCoordm)[i] /= myTopo->atoms[i].scaledMass;
    }
  }

  // fluctuation using Dr. Skeel's LI scheme which involves a semi-update
  // of velocities and a complete update of positions
  void NormalModeUtilities::nmlDrift(Vector3DBlock *myPositions, Vector3DBlock *myVelocities, Real dt, GenericTopology *myTopo) {

    //const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real tau1 = ( 1.0 - exp( -myGamma * dt ) ) / myGamma;
    const Real tau2 = ( 1.0 - exp( -2 * myGamma * dt ) ) / ( 2 * myGamma );
    const Real forceConstant = 2 * Constant::BOLTZMANN * myTemp * myGamma;
    const Real sqrtTau2 = sqrt( tau2 );

    //  It is possible that roundoff and/or truncation can make this value < 0.
    //  It should be set to 0 if this happens.
    Real sqrtVal1 = ( dt - tau1 * tau1 / tau2 );

    if( sqrtVal1 < 0. )
        sqrtVal1 = 0;
    else
        sqrtVal1 = sqrt( sqrtVal1 );

    Real sqrtFCoverM = sqrt( forceConstant );/// mass );
    Real langDriftVal = sqrtFCoverM / myGamma;
    Real langDriftZ1 = langDriftVal * ( tau1 - tau2 ) / sqrtTau2;
    Real langDriftZ2 = langDriftVal * sqrtVal1;

    //generate 1st and 2nd set of random force variables and project into sub space
    if(genCompNoise) genProjGaussC(&gaussRandCoord1, &tempV3DBlk, myTopo); //and comliment if required
    else genProjGauss(&gaussRandCoord1, myTopo);
    genProjGauss(&gaussRandCoord2, myTopo);
    // Update positions and correct (semi-update) velocity
    for( unsigned int i = 0; i < myPositions->size(); i++ ) {
      posTemp[i] =(gaussRandCoord1[i]*langDriftZ1 +gaussRandCoord2[i]*langDriftZ2 +(*myVelocities)[i])*tau1;
      // semi-update velocities
      (*myVelocities)[i] = (*myVelocities)[i]*exp(-myGamma*dt)+gaussRandCoord1[i]*sqrtFCoverM*sqrtTau2;
    }
    //Check position change is in the sub space and add to positions
    subspaceVelocity(&posTemp, &posTemp);
    for( unsigned int i = 0; i < myPositions->size(); i++ )
            (*myPositions)[i]+=posTemp[i];
    //fix COM/momentum (not conserved)
    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);
  }

  //*************************************************************************************
  //****Diagonalization/ minimizer etc.**************************************************
  //*************************************************************************************

  //Diagonalize hessian******************************************************************//
  int NormalModeUtilities::diagHessian(double *eigVecO, double *eigValO, double *hsnhessM, int dim, int &numFound){
	double *wrkSp;
	int *isuppz, *iwork;

	wrkSp = new double[26*dim];
	isuppz = new int[2*dim];
	iwork = new int[10*dim];
	//Diagonalize
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
    char cmach = 's'; //String is actualy safe min but it is shortened to remove the warning
    int info = 0;				/* output 0=success */
	int m = 0;
    //call LAPACK
    abstol = Lapack::dlamch( &cmach );	//find machine safe minimum

    Lapack::dsyevr( &jobz, &range, &uplo, &n, hsnhessM, &lda, &vl, &vu, &il, &iu, &abstol, &m, eigValO, eigVecO, &ldz, isuppz,
                wrkSp, &lwork, iwork, &liwork, &info);

	numFound = m;

    //delete arrays
    delete [] iwork;
    delete [] isuppz;
    delete [] wrkSp;

    //return status
    return info;
  }

  //sort vectors for absolute value*******************************************************//
  void NormalModeUtilities::absSort(double *eigVec, double *eigVal, int *eigIndx, int dim){
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

  //find Hessian projected into mxm 'inner' space
  void NormalModeUtilities::getInnerHess(double *eigVec, double *hess, double *innerHess){
    char transA = 'T'; char transB = 'N';
    int m = _rfM; int n = _3N; int k = _3N;
    int lda = _3N; int ldb = _3N; int ldc = _rfM;
    double alpha = 1.0;	double beta = 0.0;
    double *tempMat = new double[_rfM*_3N];

    Lapack::dgemm(&transA, &transB, &m, &n, &k, &alpha, eigVec, &lda, hess, &ldb, &beta, tempMat, &ldc);
    n = _rfM; lda = _rfM;
    Lapack::dgemm(&transB, &transB, &m, &n, &k, &alpha, tempMat, &lda, eigVec, &ldb, &beta, innerHess, &ldc);

    delete [] tempMat;
  }

  //product of eigenvectors and diagonalized 'inner' hessian
  void NormalModeUtilities::getNewEigs(double *eigVec, double *origEigVec, double *innerEigVec){
    char transA = 'N'; char transB = 'N';
    int m = _3N; int n = _rfM; int k = _rfM;
    int lda = _3N; int ldb = _rfM; int ldc = _3N;
    double alpha = 1.0;	double beta = 0.0;

    Lapack::dgemm(&transA, &transB, &m, &n, &k, &alpha, origEigVec, &lda, innerEigVec, &ldb, &beta, eigVec, &ldc);

    //copy back to eigVec
    //for(int i=0;i<_rfM*_3N;i++) eigVec[i] = tempMat[i];
    //delete [] tempMat;
  }

  //calculate the Rayleigh quotient
  double NormalModeUtilities::calcRayleigh(double *rQ, double *boundRq, double *hsnhessM, int numv, double raylAverage){

    //Norm Mode diagnostics-Rayleigh Quotient
    //Get current Hessian
    for(int i=0 ; i<_3N*_3N ; i++) hsnhessM[i] /= raylAverage;	//divide sum of Hessians
    //set BLAS variables
    char transA = 'N';
    int m = _3N; int n = _3N;
    int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;
    //calculate Rayleigh Quo/bound
    //double rQ, boundRq;
    //define temporary position/force vector
    double *tmpFX = new double[_3N];

    Lapack::dgemv(&transA, &m, &n, &alpha, hsnhessM, &m, &((*Q)[_3N*(numv)]), &incxy, &beta, tmpFX, &incxy);
    *rQ = Lapack::ddot(&n, &((*Q)[_3N*(numv)]), &incxy, tmpFX, &incxy);

    //bound
    for(int i=0;i<_3N;i++) hsnhessM[i+_3N*i] -= *rQ;

    Lapack::dgemv(&transA, &m, &n, &alpha, hsnhessM, &m, &((*Q)[_3N*(numv)]), &incxy, &beta, tmpFX, &incxy);
    *boundRq = Lapack::dnrm2(&n, tmpFX, &incxy);

    for(int i=0;i<_3N;i++) hsnhessM[i+_3N*i] += *rQ; //restore

    delete [] tmpFX;

    return *rQ;
  }

  //*************************************************************************************
  //****Minimizer************************************************************************
  //*************************************************************************************

  //Simple-Steepest decent minimizer for all modes outside subspace.
  //PK update respects 'c' subspace positions. Requires virtual force calculation function utilityCalculateForces().
  int NormalModeUtilities::minimizer(Real peLim, int numloop, bool simpM, bool reDiag, bool nonSubspace, int *forceCalc, Real *lastLambda,
                                ScalarStructure *myEnergies, Vector3DBlock *myPositions, GenericTopology *myTopo) {
    int in, itr, numLambda;
    Real oldPot, lambda, lambda1, lambdaSlp, lambdaSlp1, lastDiff;
    int rsCG;

    //initialize
    rsCG = 0; *lastLambda = 0.0; numLambda = 0; itr = 0; *forceCalc = 0; lastDiff = 5; //start at 5 kcal mol^{-1}

    utilityCalculateForces();
    (*forceCalc)++;
    //Set first value of /lambda to be 1/eigval
    lambda = 1.0 / *eigValP;	//exact solution for highest frequency mode if force mass weighted
    for(in=0;in<numloop;in++){
        itr++;
        //
        report.precision(10);
        report <<debug(6)<<"[NormalModeUtilities::minimizer] PE= "<<myEnergies->potentialEnergy()<<endr;
        //****find search direction vector posTemp
        //find forces in compliment space.
        if(nonSubspace) nonSubspaceForce(myForcesP, myForcesP);
        //sift so position move is in compliment space, mass weighted. Replaces nonSubspacePosition(myForces, myForces)
        for(int i=0;i<_N;i++) (*myForcesP)[i] /= myTopo->atoms[i].scaledMass;
        //set posTemp
        posTemp = *myForcesP;
        //find slope of original PE with /lambda=0 here.
        lambdaSlp1 = 0.0;
        for( int k = 0; k < _N; k++ ) lambdaSlp1 -= posTemp[k].dot((*myForcesP)[k]);
        //save PE at /lambda=0
        oldPot = myEnergies->potentialEnergy();
        report <<debug(7)<<"[NormalModeUtilities::minimizer] lambd= "<<lambda<<endl;
        //find force at new position pos+/lambda*posTemp
        (*myPositions).intoWeightedAdd(lambda,posTemp);
        utilityCalculateForces();
        (*forceCalc)++;
        //Full minimizer? then solve for quadratic minimum
        if(!simpM){
            //find slope of PE with /lambda here
            lambdaSlp = 0.0;
            for( int k = 0; k < _N; k++ ) lambdaSlp -= posTemp[k].dot((*myForcesP)[k]);
            //solve for minimum for quadratic fit using two PE vales and the slope with /lambda
            Real a, b, oldLambda;
            oldLambda = lambda;
            a = -((myEnergies->potentialEnergy() - oldPot) / lambda - lambdaSlp) / lambda;
            b = lambdaSlp - 2.0 * a * lambda;
            lambda = -b / (2 * a);
            if(lambda <= 0.0) lambda = oldLambda;
            else{
                //Put solution into positions (but remove temporary solution for quadratic fit via oldLambda)
                (*myPositions).intoWeightedAdd(lambda-oldLambda,posTemp);
                utilityCalculateForces();
                (*forceCalc)++;
            }
        }
        //end full minimizer
        //update total gamma
        *lastLambda += lambda;
        numLambda++;
        //test for end, too large lambda test first
        if((oldPot - myEnergies->potentialEnergy()) < 0){
          if(rsCG>4){
            report << error << "[NormalModeUtilities::minimizer] Minimization failed, Aborting. "<<rsCG <<endr;
          }else{
              if(!reDiag){  //allow minimization if mode at angle to sub-space
                  //calc optimum lambda from first slope
                  Real a1;
                  a1 = (myEnergies->potentialEnergy() - oldPot - lambdaSlp1 * lambda) / (lambda * lambda);
                  lambda1 = -lambdaSlp1 / (2 * a1);
                  //Test that the quadratic solution gives a predicted PE value
                  //where the difference from the old PE value is bounded by the difference
                  //from the last succesful step, else solve quadratic for the last difference.
                  Real calcPE = a1*lambda1*lambda1+lambdaSlp1*lambda1+oldPot;
                  if(oldPot - calcPE > lastDiff && lambdaSlp1 != 0.0){
                      lambda1 = (-lastDiff * 2.0) / lambdaSlp1;
                  }
                  (*myPositions).intoWeightedAdd(-lambda,posTemp);		//reset positions
                  *lastLambda -= lambda;
                  numLambda--;
                  utilityCalculateForces();
                  (*forceCalc)++;
                  if(lambda1 > 0.0 && lambda1 < lambda) lambda = lambda1;
                  else lambda /= 2.0;
                  rsCG++;
                  report <<debug(1)<<"[NormalModeUtilities::minimizer] Reset CG, PE fail. Cycle= "<<rsCG<<" lambda= "<<lambda<<endl;
              }else{
                  (*myPositions).intoWeightedAdd(-lambda,posTemp);		//reset positions
                  report <<debug(1)<<"[NormalModeUtilities::minimizer] REDIAGONALIZING! PE fail  lambda= "<<lambda<<endl;
                  return -1;				//flag aborted

              }
            }
        }else{
            rsCG = 0;
            lambda = 1.0 / *eigValP;	//revert to original value
        }
        if((oldPot - myEnergies->potentialEnergy()) < peLim && !rsCG) break;
        if(!rsCG) lastDiff = oldPot - myEnergies->potentialEnergy();
        //
    }
    if(numLambda) *lastLambda /= (Real)numLambda;
    else lastLambda = 0;
    return itr;
  }


}
