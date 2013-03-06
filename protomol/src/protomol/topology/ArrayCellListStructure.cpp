#include <protomol/topology/ArrayCellListStructure.h>
#include <protomol/base/Exception.h>
#include <protomol/base/Report.h>
#include <exception>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

ArrayCellListStructure::ArrayCellListStructure() :
	valid( false ), myArray( ArraySizes( 0 )( 0 )( 0 ) ),
	myCellSize( Vector3D( 0.0, 0.0, 0.0 ) ), myNX( 0 ), myNY( 0 ), myNZ( 0 ),
	myBegin( myArray.begin() ), myEnd( myArray.end() ),
	myBeginConst( myArray.begin() ), myEndConst( myArray.end() ), mySize( 0 ), bInit( false ), mVolume( 0.0f ),
	mMaxVolume( 16.0f ){
}

void ArrayCellListStructure::MaximumVolume( const float value ){
	mMaxVolume = value;
}

void ArrayCellListStructure::initialize( const Vector3D &max, Vector3D cellSize ) {
	if( !bInit ) {
		bInit = true;
		mVolume = max[0] * max[1] * max[2];
	} else {
		float vol = max[0] * max[1] * max[2];
		float volDiff = std::fabs( vol / mVolume );

#ifdef DEBUG_CELLVOLUME
		std::cout << "Original: " << mVolume << " New: " << vol << " Diff: " << volDiff << std::endl;
#endif

		const float maxVolumeInverse = 1.0f / (float)mMaxVolume;
		if( volDiff > mMaxVolume || volDiff < maxVolumeInverse ) {
			report << error << "Simulation volume shrunk or grew by over " << mMaxVolume << "x. Assuming simulation failure." << endr;
		}
	}

	int nx = std::max( 1, ( int )floor( max.c[0] / cellSize.c[0] + Constant::EPSILON ) );
	int ny = std::max( 1, ( int )floor( max.c[1] / cellSize.c[1] + Constant::EPSILON ) );
	int nz = std::max( 1, ( int )floor( max.c[2] / cellSize.c[2] + Constant::EPSILON ) );

	if( nx != myNX || ny != myNY || nz != myNZ || myCellSize != cellSize ) {
#ifdef DEBUG_ARRAYCELLLISTSTRUCTURE
		report << hint << "CubicCellManager: re-size from (" << myNX << "," <<
			   myNY << "," << myNZ << ") to (" << nx << "," << ny << "," << nz <<
			   ")." <<
			   endr;
		report << hint << "CubicCellManager: N_x " <<
			   toString( max.x / cellSize.x ) << "," << ( int )floor(
				   max.x / cellSize.x + Constant::EPSILON ) << "," << ( int )floor(
				   max.x / cellSize.x + 1 - Constant::EPSILON ) << "," << ( int )floor(
				   max.x / cellSize.x + 1 + Constant::EPSILON ) << endr;
#endif

		if( nx * ny * nz * sizeof( T ) >= power < sizeof( size_t ) * 8 > ( 2.0 ) )
			report << error
				   << "Your systems is expanding such that the "
				   << "the ratio simulation box / cell size is to big. "
				   << "You may decrease your timestep in your integrator or "
				   << "increase your cell size. End of advice."
				   << endr;
		myCellSize = cellSize;
		myNX = nx;
		myNY = ny;
		myNZ = nz;

		try {
			myArray.resize( ArraySizes( myNX )( myNY )( myNZ ) );
		} catch( std::bad_alloc &e ) {
			report << error <<
				   "[ArrayCellListStructure::initialize] Could not allocate memory for "
				   "CellListStructure[" << myNX << "][" << myNY << "][" << myNZ << "]."
				   << '\n' << "Try adjusting your integration parameters (time step, "
				   << "number of modes (NML), etc.)" << endr;
		}

		int max2 = Constant::MAX_INT_2;
		myMaxNX = max2 - max2 % nx;
		myMaxNY = max2 - max2 % ny;
		myMaxNZ = max2 - max2 % nz;
		myNX1 = -myNX / 2;
		myMaxNX1 = myMaxNX - myNX1;
		myNY1 = -myNY / 2;
		myMaxNY1 = myMaxNY - myNY1;
		myNZ1 = -myNZ / 2;
		myMaxNZ1 = myMaxNZ - myNZ1;

		for( int i = 0; i < myNX; i++ ) {
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
			Array<T, 3>::RefArray<2> z2 = myArray[i];
#else
			RefArray<T, 2> z2 = myArray[i];
#endif
			for( int j = 0; j < myNY; j++ ) {
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
				Array<T, 3>::RefArray<1> z1 = z2[j];
#else
				RefArray<T, 1> z1 = z2[j];
#endif
				for( int k = 0; k < myNZ; k++ ) {
					z1[k].first = T1( i, j, k );
					z1[k].second = -1;
				}
			}
		}
	} else {
		T *a = myArray.begin();
		for( unsigned int i = 0; i < myArray.size(); i++ ) {
			a[i].second = -1;
		}
	}

	valid = false;
	myBegin = myArray.begin();
	myEnd = myArray.end();
	myBeginConst = myArray.begin();
	myEndConst = myArray.end();
	mySize = myArray.size();

#ifdef DEBUG_ARRAYCELLLISTSTRUCTURE
	report << plain
		   << "[ArrayCellListStructure::initialize] cellsize=" << cellSize <<
		   ", box(" << myNX << "," << myNY << "," << myNZ << ")."
		   << endr;
#endif
}

void ArrayCellListStructure::updateCache() {
	T *a = myArray.begin();
	mySize = 0;
	myBegin = myArray.begin();
	myEnd = myArray.begin();
	myBeginConst = myArray.begin();
	myEndConst = myArray.begin();
	unsigned int count = myArray.size();
	bool first = true;
	for( unsigned int i = 0; i < count; i++ )
		if( a[i].second >= 0 ) {
			if( first ) {
				myBegin = myArray.begin() + i;
				myBeginConst = myArray.begin() + i;
			}
			myEnd = myArray.begin() + i + 1;
			myEndConst = myArray.begin() + i + 1;
			first = false;
			mySize++;
		}

	valid = true;
}

