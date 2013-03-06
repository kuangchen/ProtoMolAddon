#ifndef CUBICCELLMANAGER_H
#define CUBICCELLMANAGER_H

#include <protomol/topology/ArrayCellListStructure.h>
#include <protomol/config/Parameter.h>
#include <protomol/type/Vector3D.h>

#include <vector>

namespace ProtoMol {
	/**
	 * The cell manager for equal-sized (cubic) cells. For optimization reasons
	 * in case of periodic boundary conditions the cells are not cubic any more
	 * in order to fit the system by multiples of the cell dimensions.
	 */
	class CubicCellManager {
		public:
			/// topology and cell location structure of the cell
			typedef CubicCellLocation Cell;
			/// implementation of the cell list
			typedef ArrayCellListStructure CellListStructure;
		public:
			static const std::string keyword;
		private:
			Real mCellSize;
			mutable Vector3D mCellSizeVector;
			mutable Vector3D mCellSizeVectorInverse;
		public:
			CubicCellManager();
			CubicCellManager( const Real cellsize );
		public:
			Real getCellSize() const;
			Real getCellVolume() const;
			Vector3D getCellSizeVector() const;

			Cell findCell( const Vector3D &position ) const;

			void initialize( CellListStructure &cellList, const Vector3D &min, const Vector3D &max, const bool pbc ) const;
			void updateCache( CellListStructure &cellList ) const;

			const std::string &getKeyword() const;
			void getParameters( std::vector<Parameter> &parameters ) const;
			static CubicCellManager make( const std::vector<Value>& values );
	};
}
#endif /* CUBICCELLMANAGER_H */
