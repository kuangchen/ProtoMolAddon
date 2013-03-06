#include "protomol/type/EigenvectorInfo.h"

using namespace ProtoMol;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructors, destructors, assignment
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EigenvectorInfo::EigenvectorInfo() : myEigenvectorLength( 0 ), myNumEigenvectors( 0 ),
	myNumUsedEigenvectors( 0 ),
	myEigenvectors( 0 ), myOrigCEigval( 0.0 ), myNewCEigval( 0.0 ), myOrigTimestep( 0.0 ),
	reDiagonalize( false ), havePositionsChanged( false ), OpenMMMinimize( false ), RediagonalizationCount( 0 ),
	mySingleEigs( 0 ), myEigVecChanged( true ), myMinimumLimit( 0.5 ), currentMode( -1 ) {

}

EigenvectorInfo::EigenvectorInfo( unsigned int n, unsigned int m ) : myEigenvectorLength( n ),
	myNumEigenvectors( m ), myNumUsedEigenvectors( 0 ),
	myMaxEigenvalue( 0.0 ), myEigenvectors( new double[n *m * 3] ),
	myOrigCEigval( 0.0 ), myNewCEigval( 0.0 ), myOrigTimestep( 0.0 ), reDiagonalize( false ),
	havePositionsChanged( false ), OpenMMMinimize( false ), RediagonalizationCount( 0 ), mySingleEigs( 0 ), myEigVecChanged( true ), myMinimumLimit( 0.5 ),
	currentMode( -1 ) {}

EigenvectorInfo::~EigenvectorInfo() {
	if( myEigenvectors ) {
		delete [] myEigenvectors;
	}
	if( mySingleEigs ) {
		delete [] mySingleEigs;
	}
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New methods of class EigenvectorInfo
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool EigenvectorInfo::initializeEigenvectors() {
	try {
		myEigenvectors = new double[myEigenvectorLength * myNumEigenvectors * 3];
	} catch( std::bad_alloc & ) {
		return false;
	}

	return true;
}

float* EigenvectorInfo::getFloatEigPointer() {
	//no eigenvectors assigned?
	if( !myEigenvectors ) {
		return ( float * )0;
	}

	const unsigned int arrayLen = myEigenvectorLength * myNumEigenvectors * 3;

	//assign storage if required
	if( !mySingleEigs ) {
		try {
			mySingleEigs = new float[arrayLen];
		} catch( std::bad_alloc & ) {
			return ( float * )0;
		}
	}

	//update values if double array updated
	if( myEigVecChanged ) {
		for( unsigned int i = 0; i < arrayLen; i++ ) {
			mySingleEigs[i] = ( float )myEigenvectors[i];
		}

		myEigVecChanged = false;
	}

	return mySingleEigs;
}
