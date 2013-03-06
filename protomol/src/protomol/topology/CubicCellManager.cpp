#include <protomol/topology/CubicCellManager.h>

#include <protomol/base/Report.h>
#include <protomol/base/Exception.h>

namespace ProtoMol {
	using namespace ProtoMol::Report;

	/* Local Helper Functions */
	const Vector3D InverseOf( const Vector3D &vector ) {
		return Vector3D( 1.0 / vector.c[0], 1.0 / vector.c[1], 1.0 / vector.c[2] );
	}

	/* Class Static Variables */
	const string CubicCellManager::keyword( "Cubic" );

	/* Class Functions */
	CubicCellManager::CubicCellManager() : mCellSize( 0.0 ) {

	}

	CubicCellManager::CubicCellManager( const Real cellSize )
		: mCellSize( cellSize ), mCellSizeVector( Vector3D( cellSize, cellSize, cellSize ) ),
		  mCellSizeVectorInverse( Vector3D( 1.0 / cellSize, 1.0 / cellSize, 1.0 / cellSize ) ) {

	}

	Real CubicCellManager::getCellSize() const {
		return mCellSize;
	}

	Real CubicCellManager::getCellVolume() const {
		return mCellSizeVector.c[0] * mCellSizeVector.c[1] * mCellSizeVector.c[2];
	}

	Vector3D CubicCellManager::getCellSizeVector() const {
		return mCellSizeVector;
	}

	CubicCellManager::Cell CubicCellManager::findCell( const Vector3D &position ) const {
		return Cell( ( int )floor( position.c[0] * mCellSizeVectorInverse.c[0] ),
					 ( int )floor( position.c[1] * mCellSizeVectorInverse.c[1] ),
					 ( int )floor( position.c[2] * mCellSizeVectorInverse.c[2] ) );
	}

	void CubicCellManager::initialize( CellListStructure &cellList,
									   const Vector3D &min, const Vector3D &max,
									   const bool pbc ) const {

		const Vector3D vDiff( max - min );
		Vector3D vResult( vDiff );

		if( pbc ) {
			const Real x = vDiff.c[0] / std::max( 1.0, floor( vDiff.c[0] / mCellSize ) );
			const Real y = vDiff.c[1] / std::max( 1.0, floor( vDiff.c[1] / mCellSize ) );
			const Real z = vDiff.c[2] / std::max( 1.0, floor( vDiff.c[2] / mCellSize ) );

			mCellSizeVector = Vector3D( x, y, z );
			mCellSizeVectorInverse = InverseOf( mCellSizeVector );
		} else {
			vResult = vResult + mCellSizeVector;
		}

		cellList.initialize( vResult, mCellSizeVector );
	}

	void CubicCellManager::updateCache( CellListStructure &cellList ) const {
		cellList.updateCache();
	}

	const std::string &CubicCellManager::getKeyword() const {
		return keyword;
	}

	void CubicCellManager::getParameters( std::vector<Parameter> &parameters ) const {
		parameters.push_back
		( Parameter( "cellSize", Value( mCellSize, ConstraintValueType::Positive() ),
					 Text( "For Periodic BC this must be < least cell basis vector."
						   "  Typically 1/2 of the least cutoff value" ) ) );
		parameters.push_back
		( Parameter( "maxvolume", Value( 16, ConstraintValueType::Positive() ),
					 Text( "Maximum volume change for cell" ) ) );
	}

	CubicCellManager CubicCellManager::make( const std::vector<Value>& values ) {
		Real cellSize = 1.0;

		if( !values[0].get( cellSize ) || cellSize <= 0.0 ) {
			report	<< error << keyword
					<< " cellmanager: cutoff > 0 (" << values[0].getString() << ")." << endr;
		}

		return CubicCellManager( cellSize );
	}
}
