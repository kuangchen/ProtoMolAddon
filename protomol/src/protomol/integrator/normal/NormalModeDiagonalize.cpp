#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/BlockMatrix.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/STSIntegrator.h>

#include<iostream>
#include<fstream>

#include <protomol/base/Lapack.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

const string NormalModeDiagonalize::keyword( "NormalModeDiagonalize" );

NormalModeDiagonalize::NormalModeDiagonalize() :
	MTSIntegrator(), NormalModeUtilities(), firstDiag( true ),
	fullDiag( 0 ), removeRand( 0 ), rediagCount( 0 ), nextRediag( 0 ),
	validMaxEigv( 0 ), myNextNormalMode( 0 ), myLastNormalMode( 0 ),
	rediagHysteresis( 0 ), hessianCounter( 0 ), rediagCounter( 0 ),
	rediagUpdateCounter( 0 ), eigenValueThresh( 0 ), blockCutoffDistance( 0 ),
	blockVectorCols( 0 ), residuesPerBlock( 0 ), memory_Hessian( 0 ),
	memory_eigenvector( 0 ), checkpointUpdate( false ), origCEigVal( 0 ),
	origTimestep( 0 ), autoParmeters( false ), adaptiveTimestep( 0 ),
	postDiagonalizeMinimize( 0 ), minLim( 0 ), maxMinSteps( 0 ),
	geometricfdof( false ), numerichessians( false ),
	mShouldDoMultipleRediagonalization( false ), mMaxRediagonalizations( 5 )  {
}

NormalModeDiagonalize::
NormalModeDiagonalize( int cycles, int redi, bool fDiag, bool rRand,
					   Real redhy, Real eTh, int bvc, int rpb, Real dTh,
					   bool apar, bool adts, bool pdm, Real ml, int maxit,
					   bool geo, bool num, bool multRediag, unsigned int maxRediags,
					   ForceGroup *overloadedForces,
					   StandardIntegrator *nextIntegrator ) :
	MTSIntegrator( cycles, overloadedForces, nextIntegrator ),
	NormalModeUtilities( 1, 1, 91.0, 1234, 300.0 ), firstDiag( true ),
	fullDiag( fDiag ), removeRand( rRand ), rediagCount( redi ), nextRediag( 0 ),
	validMaxEigv( 0 ), myNextNormalMode( 0 ), myLastNormalMode( 0 ),
	rediagHysteresis( redhy ), hessianCounter( 0 ), rediagCounter( 0 ),
	rediagUpdateCounter( 0 ), eigenValueThresh( eTh ),
	blockCutoffDistance( dTh ), blockVectorCols( bvc ),
	residuesPerBlock( rpb ),  memory_Hessian( 0 ), memory_eigenvector( 0 ),
	checkpointUpdate( false ), origCEigVal( 0 ), origTimestep( 0 ),
	autoParmeters( apar ), adaptiveTimestep( adts ), postDiagonalizeMinimize( pdm ),
	minLim( ml ), maxMinSteps( maxit ), geometricfdof( geo ), numerichessians( num ),
	mShouldDoMultipleRediagonalization( multRediag ), mMaxRediagonalizations( maxRediags ) {

	//find forces and parameters
	rHsn.findForces( overloadedForces );
}

NormalModeDiagonalize::~NormalModeDiagonalize() {
	//output stats
	report.precision( 5 );
	if( rediagCounter && hessianCounter ) {
		report << plain
			   << "NML Timing: Hessian: " << ( blockDiag.hessianTime.getTime() ).getRealTime() << "[s] (" << hessianCounter << " times)"
			   << " diagonalize: " << ( blockDiag.rediagTime.getTime() ).getRealTime() << "[s] (" << rediagCounter << " re-diagonalizations)."
			   << endl;

		if( !fullDiag ) {
			report << plain << "NML Memory: "
				   << "Hessian: "     << memory_Hessian     << "[Mb], "
				   << "diagonalize: " << memory_eigenvector << "[Mb], "
				   << "vectors: "     << _3N *_rfM *sizeof( double ) / 1000000 << "[Mb]."
				   << endl;
		}
	}
}

void NormalModeDiagonalize::initialize( ProtoMolApp *app ) {
	MTSIntegrator::initialize( app );

	// Setup adaptive timestep data from checkpoint
	app->eigenInfo.myOrigCEigval = origCEigVal;
	app->eigenInfo.myOrigTimestep = origTimestep;

	//test topmost integrator
	if( top() != this ) {
		report << error << "NormalModeDiagonalize not top integrator." << endr;
	}

	myNextNormalMode  = dynamic_cast<NormalModeUtilities *>( myNextIntegrator );

	//Using complement of next integrator, so copy
	firstMode = myNextNormalMode->firstMode;
	numMode   = myNextNormalMode->numMode;

	//NM initialization, but fix dof, as set by next integrator, use complementry forces and dont generate noise
	NormalModeUtilities::initialize( ( int )app->positions.size(), app, myForces, COMPLIMENT_FORCES );
	_rfM = myNextNormalMode->_rfM;

	myLastNormalMode  = dynamic_cast<NormalModeUtilities *>( bottom() );

	//do first force calculation, and remove non sub-space part
	app->energies.clear(); //Need this or initial error, due to inner integrator energy?
	initializeForces();

	//Initialize Hessian array, OR assign hessian array for residues.
	if( fullDiag ) {
		rHsn.initialData( _3N );
	} else {
		//automatically generate parameters?
		if( autoParmeters ) {
			residuesPerBlock = ( int )pow( ( double )_N, 0.6 ) / 15;
			blockVectorCols = 10 + ( int )sqrt( ( float )residuesPerBlock );
			blockCutoffDistance = rHsn.cutOff;
			report << debug( 1 ) << "[NormalModeDiagonalize::initialize] Auto parameters: residuesPerBlock " << residuesPerBlock <<
				   ", blockVectorCols " << blockVectorCols <<
				   ", blockCutoffDistance " << blockCutoffDistance << "." << endr;
		}

		//assign hessian array for residues, and clear.
		rHsn.initialResidueData( app->topology, residuesPerBlock, ( blockCutoffDistance == 0.0 ) );
	}

	// Check if array is already assigned
	if( app->eigenInfo.myEigenvectors ) {
		firstDiag = false;
		validMaxEigv = true;
	} else {
		firstDiag = true;
		validMaxEigv = false;

		//Calculate array size to be created
		app->eigenInfo.myEigenvectorLength = _N;
		app->eigenInfo.myNumEigenvectors = ( fullDiag == true ) ? _3N : _rfM;
		if( !app->eigenInfo.initializeEigenvectors() ) {
			report << error << "Eigenvector array allocation error." << endr;
		}

	}

	//flag used eigs
	app->eigenInfo.myNumUsedEigenvectors = _rfM;

	//Initialize BlockHessianDiagonalize, pass BlockHessian if Blocks (not full diag)
	if( fullDiag ) {
		blockDiag.initialize( _3N );
	} else {
		blockDiag.initialize( &rHsn, _3N, ( StandardIntegrator * )this );
	}

	//Diagnostics
	memory_Hessian = memory_eigenvector = 0;

	//setup rediag counter in case valid
	nextRediag = app->currentStep;//( int )( app->topology->time / getTimestep() ); //rediag first time

	//save positions where diagonalized for checkpoint save (assume I.C. if file)
	if( !checkpointUpdate ) {
		diagAt = app->positions;
	} else {
		firstDiag = true;
	}

	newDiag = false;

	//timers/counters for diagnostics
	hessianCounter = rediagCounter = rediagUpdateCounter = 0;
}

//*************************************************************************************
//****Normal run routine***************************************************************
//*************************************************************************************

void NormalModeDiagonalize::run( int numTimesteps ) {
	if( numTimesteps < 1 ) {
		return;
	}

	//Current step at start
	int currentStepNum = app->currentStep;

	//main loop
	app->energies.clear();
	for( int i = 0; i < numTimesteps; ) {

		//Diagonalization if repetitive, first for forced
		if( ( rediagCount && currentStepNum >= nextRediag ) || firstDiag || app->eigenInfo.reDiagonalize ) {

			if( ( rediagCount && currentStepNum >= nextRediag ) || firstDiag ){
				nextRediag += rediagCount;
			}

			newDiag = true;

			report << debug( 2 ) << "[NormalModeDiagonalize::run] Finding diagonalized Hessian." << endr;

			if( !( checkpointUpdate && firstDiag ) ) {
				//remove last random perturbation?
				if( removeRand ) {
					app->positions.intoSubtract( myLastNormalMode->gaussRandCoord1 );
				}
			}

			//Diagonalize
			if( fullDiag ) {
				FullDiagonalize();
			} else {
				CoarseDiagonalize();
			}
			
			if( !( checkpointUpdate && firstDiag ) ) {
				diagAt = app->positions;
				//remove last random perturbation?
				if( removeRand && !postDiagonalizeMinimize && fullDiag ) {
					app->positions.intoAdd( myLastNormalMode->gaussRandCoord1 );
				}
			}
			
			//adaptive timestep?
			const double newCEig = app->eigenInfo.myNewCEigval;
			const double oldCEig = app->eigenInfo.myOrigCEigval;
			const double baseTimestep = app->eigenInfo.myOrigTimestep;
			
			if( adaptiveTimestep && baseTimestep > 0 && newCEig > 0 && newCEig != oldCEig ) {
				
				const double tRatio = sqrt( oldCEig / newCEig );
				
				const double oldTimestep = bottom()->getTimestep();
				
				if( baseTimestep *tRatio <= baseTimestep ) {
					( ( STSIntegrator * )bottom() )->setTimestep( baseTimestep * tRatio );
					
					report << debug( 1 ) << "Adaptive time-step change, base " << baseTimestep <<
					", new " << baseTimestep *tRatio << ", old " << oldTimestep << "." << endr;
				}
				
			}
			
			//sift current velocities/forces
			myNextNormalMode->subSpaceSift( &app->velocities, myForces );
			
			//clear re-diag flag
			app->eigenInfo.reDiagonalize = false;
		}

		//run integrator
		int stepsToRun = min( nextRediag - currentStepNum, numTimesteps - i );

		myNextIntegrator->run( stepsToRun );
		i += stepsToRun;
		currentStepNum += stepsToRun;

		//remove diagonalization flags after inner integrator call
		newDiag = false;
	}

}

bool NormalModeDiagonalize::Minimize(){
	if( !app->eigenInfo.OpenMMMinimize ) {
		Real lastLambda; int forceCalc = 0; //diagnostic/effective gamma
		
		//do minimization with local forces, max loop maxMinSteps, set subSpace minimization true
		int itrs = minimizer( minLim, maxMinSteps, true, false, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology );
		app->eigenInfo.havePositionsChanged = true;
		
		report << debug( 2 ) << "[NormalModeDiagonalize::run] iterations = " << itrs << " force calcs = " << forceCalc << endr;
		
		// Break if termination condition is met
		if( itrs <= 2 ) {
			return true;
		}
	}
	
	return true;
}

void NormalModeDiagonalize::FullDiagonalize() {
	//****Full method**********************************************************************//
	// Uses BLAS/LAPACK to do 'brute force' diagonalization                                //
	//*************************************************************************************//
	report << debug( 2 ) << "Start diagonalization." << endr;
	
	int loops = 1;
	if( app->eigenInfo.reDiagonalize || mShouldDoMultipleRediagonalization ) {
		loops = mMaxRediagonalizations;
	}
	for( int iteration = 0; iteration < loops; iteration++ ) {
		//Find Hessians
		blockDiag.hessianTime.start(); //time Hessian
		rHsn.clear();
		rHsn.evaluate( &app->positions, app->topology, true ); //mass re-weighted hessian, used to pass diagAt, but this defeats the multiple diagonalize
		report << debug( 2 ) << "Hessian found." << endr;
		
		//stop timer
		blockDiag.hessianTime.stop();
		hessianCounter++;
		
		//Diagonalize
		blockDiag.rediagTime.start();
		int numeFound;
		int info = blockDiag.diagHessian( *Q , blockDiag.eigVal, rHsn.hessM, _3N, numeFound );
		
		if( info ) {
			report << error << "Full diagonalization failed." << endr;
		}
		
		//find number of -ve eigs
		int ii;
		for( ii = 0; ii < _3N - 3; ii++ ) {
			if( blockDiag.eigVal[ii + 3] > 0 ) {
				break;
			}
		}
		
		report << debug( 1 ) << "[NormalModeDiagonalize::run] Full diagonalize. No. negative eigenvales = " << ii << endr;
		
		for( int i = 0; i < _3N; i++ ) {
			blockDiag.eigIndx[i] = i;
		}
		
		blockDiag.absSort( *Q , blockDiag.eigVal, blockDiag.eigIndx, _3N );
		
		//set new max eigenvalue in C
		app->eigenInfo.myNewCEigval = fabs( blockDiag.eigVal[_rfM] ); //safe as eigval set t length sz=_3N >= _rfM
		
		//flag update to eigenvectors
		*eigVecChangedP = true;
		
		blockDiag.rediagTime.stop();
		rediagCounter++;
		
		//set flags if firstDiag
		if( firstDiag ) {
			numEigvectsu = _3N;
			*eigValP = blockDiag.eigVal[_3N - 1];
			
			//first max eigenvalue in C, save original timestep for adaptive use
			if( !checkpointUpdate ) {
				app->eigenInfo.myOrigCEigval = app->eigenInfo.myNewCEigval;
				app->eigenInfo.myOrigTimestep = bottom()->getTimestep();
			}
			
			validMaxEigv = true;
			firstDiag = false;
		}
		
		//post diag minimize?
		if( loops > 1 || postDiagonalizeMinimize ) {
			if( Minimize() ) break;
		}
	}
}

void NormalModeDiagonalize::CoarseDiagonalize(){
	//****Coarse method**************************************************************************//
	// Process:  Finds isolated 'minimized' block (of residues) Hessians   [evaluateResidues]    //
	//           Diagonalizes blocks to form block eigenvectors B          [findCoarseBlockEigs] //
	//           Finds actual Hessian H (but coarse grained) then S=B^THB  [innerHessian]        //
	//           Diagonalizes S to get eigenvectors Q, then approximate                          //
	//           eigenvectors are the first 'm' columns of BQ.                                   //
	//*******************************************************************************************//
	report << debug( 2 ) << "Start coarse diagonalization." << endr;
	
	int loops = 1;
	if( app->eigenInfo.reDiagonalize || mShouldDoMultipleRediagonalization ) {
		loops = mMaxRediagonalizations;
	}
	for( int iteration = 0; iteration < loops; iteration++ ) {
		Real max_eigenvalue = blockDiag.findEigenvectors( &app->positions, app->topology,
														 *Q , _3N, _rfM,
														 blockCutoffDistance, eigenValueThresh, blockVectorCols,
														 geometricfdof, numerichessians );
		
		//Stats/diagnostics
		rediagCounter++; hessianCounter++;
		memory_Hessian = ( rHsn.memory_base + rHsn.memory_blocks ) * sizeof( Real ) / 1000000;
		memory_eigenvector = blockDiag.memory_footprint * sizeof( Real ) / 1000000;
		
		//set new max eigenvalue in C
		app->eigenInfo.myNewCEigval = fabs( blockDiag.eigVal[_rfM] ); //safe as eigval set t length sz=_3N >= _rfM
		
		//flag update to eigenvectors
		*eigVecChangedP = true;
		
		//set flags if firstDiag (firstDiag can now be coarse)
		if( firstDiag ) {
			//Number of eigenvectors in set, _rfM
			numEigvectsu = _rfM;
			
			//use 1000 for regression tests, set REGRESSION_T NE 0.
			if( REGRESSION_T ) {
				*eigValP = 1000;
			} else {
				
				//maximum from blocks
				*eigValP = max_eigenvalue;
			}
			
			//first max eigenvalue in C, save original timestep for adaptive use
			if( !checkpointUpdate ) {
				app->eigenInfo.myOrigCEigval = app->eigenInfo.myNewCEigval;
				app->eigenInfo.myOrigTimestep = bottom()->getTimestep();
			}
			
			//flags
			validMaxEigv = true;
			firstDiag = false;
		}
		
		report << debug( 2 ) << "Coarse diagonalization complete. Maximum eigenvalue = " << max_eigenvalue << "." << endr;
		
		//post diag minimize?
		if( loops > 1 || postDiagonalizeMinimize ) {
			if( Minimize() ) break;
		}
	}
}

//********************************************************************************************************************************************

//*************************************************************************************
//****Output int paramiters************************************************************
//*************************************************************************************

void NormalModeDiagonalize::getParameters( vector<Parameter>& parameters ) const {
	MTSIntegrator::getParameters( parameters );

	parameters.push_back( Parameter( "reDiagFrequency",
									 Value( rediagCount, ConstraintValueType::NotNegative() ),
									 0,
									 Text( "Frequency of re-diagonalization (steps)." ) ) );

	parameters.push_back( Parameter( "fullDiag",
									 Value( fullDiag,            ConstraintValueType::NoConstraints() ),
									 false,
									 Text( "Full diagonalization?" ) ) );

	parameters.push_back( Parameter( "removeRand",
									 Value( removeRand,          ConstraintValueType::NoConstraints() ),
									 false,
									 Text( "Remove last random perturbation?" ) ) );

	parameters.push_back( Parameter( "rediagHysteresis",
									 Value( rediagHysteresis,    ConstraintValueType::NotNegative() ),
									 0.0,
									 Text( "Re-diagonalization hysteresis." ) ) );

	parameters.push_back( Parameter( "eigenValueThresh",
									 Value( eigenValueThresh,    ConstraintValueType::NotNegative() ),
									 5.0,
									 Text( "'Inner' eigenvalue inclusion threshold." ) ) );

	parameters.push_back( Parameter( "blockVectorCols",
									 Value( blockVectorCols,     ConstraintValueType::NotNegative() ),
									 0,
									 Text( "Target number of block eigenvector columns." ) ) );

	parameters.push_back( Parameter( "residuesPerBlock",
									 Value( residuesPerBlock,    ConstraintValueType::NotNegative() ),
									 1,
									 Text( "Residues per block." ) ) );

	parameters.push_back( Parameter( "blockCutoffDistance",
									 Value( blockCutoffDistance, ConstraintValueType::NotNegative() ),
									 10,
									 Text( "Block cutoff distance for electrostatic forces." ) ) );

	parameters.push_back( Parameter( "autoParameters",
									 Value( autoParmeters, ConstraintValueType::NoConstraints() ),
									 false,
									 Text( "Automatically generate diagonalization parameters." ) ) );

	parameters.push_back( Parameter( "adaptiveTimestep",
									 Value( adaptiveTimestep, ConstraintValueType::NoConstraints() ),
									 false,
									 Text( "Adapt time-step to latest diagonalization eigenvalues." ) ) );

	parameters.push_back( Parameter( "postDiagonalizeMinimize",
									 Value( postDiagonalizeMinimize, ConstraintValueType::NoConstraints() ),
									 false,
									 Text( "Minimize after diagonalization." ) ) );

	parameters.push_back( Parameter( "minimlim",
									 Value( minLim, ConstraintValueType::NotNegative() ),
									 0.1,
									 Text( "Minimizer target PE difference kcal mole^{-1}" ) ) );

	parameters.push_back( Parameter( "maxminsteps",
									 Value( maxMinSteps, ConstraintValueType::Positive() ),
									 100,
									 Text( "maximum number of minimizer steps." ) ) );

	parameters.push_back( Parameter( "geometricfdof",
									 Value( geometricfdof, ConstraintValueType::NoConstraints() ),
									 false, Text( "Calculate fixed degrees of freedom geometrically." ) ) );

	parameters.push_back( Parameter( "numericHessians",
									 Value( numerichessians, ConstraintValueType::NoConstraints() ),
									 false, Text( "Calculate Hessians numerically." ) ) );

	parameters.push_back( Parameter( "multiplerediag",
									 Value( mShouldDoMultipleRediagonalization, ConstraintValueType::NoConstraints() ),
									 false, Text( "Multiple Rediagonalizations?." ) ) );

	parameters.push_back( Parameter( "maxrediags",
									 Value( mMaxRediagonalizations, ConstraintValueType::NoConstraints() ),
									 5, Text( "Maximum number of rediagonalizations to do." ) ) );

}


MTSIntegrator *NormalModeDiagonalize::doMake( const vector<Value>& values, ForceGroup *fg, StandardIntegrator *nextIntegrator ) const {
	std::cout << "Value Count: " << values.size() << std::endl;
	return new NormalModeDiagonalize( values[0], values[1], values[2],
									  values[3], values[4], values[5],
									  values[6], values[7], values[8],
									  values[9], values[10], values[11],
									  values[12], values[13], values[14],
									  values[15], values[16], values[17],
									  fg, nextIntegrator );
}

//*************************************************************************************
//****Minimizers virtual force calculation*********************************************
//*************************************************************************************
void NormalModeDiagonalize::utilityCalculateForces() {
	app->energies.clear();

	calculateForces();
}

//*************************************************************************************
//****Checkpointing********************************************************************
//*************************************************************************************
void NormalModeDiagonalize::streamRead( std::istream &inStream ) {
	inStream >> diagAt;

	//adaptive timestep
	inStream >> origCEigVal;
	inStream >> origTimestep;

	checkpointUpdate = true;
}

void NormalModeDiagonalize::streamWrite( std::ostream &outStream ) const {
	outStream.precision( 15 );
	outStream << diagAt;

	//adaptive timestep
	outStream << std::endl << app->eigenInfo.myOrigCEigval;
	outStream << " " << app->eigenInfo.myOrigTimestep;
}
