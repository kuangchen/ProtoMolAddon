#ifndef ARRAYCELLSTRUCTURE_H
#define ARRAYCELLSTRUCTURE_H

#include <protomol/base/StringUtilities.h>
#include <protomol/topology/CubicCellLocation.h>
#include <protomol/type/Vector3D.h>
#include <protomol/type/Array.h>

namespace ProtoMol {
	/**
	 * Implements a cell list with an array, providing STL alike iterators and
	 * access
	 */
	class ArrayCellListStructure {
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My typedef's
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		private:
			typedef CubicCellLocation T1;
			typedef int T2;
			typedef std::pair<T1, T2> T;
			typedef Array<T, 3> TContainer;

		public:
			// STL-like types
			typedef TContainer::value_type value_type;
			typedef TContainer::reference reference;
			typedef TContainer::const_reference const_reference;
			typedef TContainer::pointer pointer;
			typedef TContainer::const_pointer const_pointer;
			typedef TContainer::iterator iterator;
			typedef TContainer::const_iterator const_iterator;
			typedef TContainer::difference_type difference_type;
			typedef TContainer::size_type size_type;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Constructors, destructors, assignment
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			ArrayCellListStructure();

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Return STL-like iterators
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			iterator begin()       {
				return myBegin;
			}
			iterator end()         {
				return myEnd;
			}
			const_iterator begin() const {
				return myBeginConst;
			}
			const_iterator end()   const {
				return myEndConst;
			}
			size_type size()  const {
				return mySize;
			}

			iterator find( const T1 &t1 ) {
				if( !checkBoundaries( t1 ) ) {
					return myEnd;
				}
				iterator itr = &myArray[t1.x][t1.y][t1.z];
				if( itr->second < 0 ) {
					return myEnd;
				}
				return itr;
			}

			const_iterator find( const T1 &t1 ) const {
				if( !checkBoundaries( t1 ) ) {
					return myEndConst;
				}
				const_iterator itr = &myArray[t1.x][t1.y][t1.z];
				if( itr->second < 0 ) {
					return myEndConst;
				}
				return itr;
			}

			T2 &operator[]( const T1 &t1 ) {
				return myArray[t1.x][t1.y][t1.z].second;
			}

			bool checkBoundaries( const T1 &t1 ) const {
				return !( t1.x < 0 || t1.y < 0 || t1.z < 0 ||
						  t1.x >= myNX || t1.y >= myNY || t1.z >= myNZ );
			}

			/// case for periodic boundary conditions
			const_iterator findPeriodic( const T1 &t1 ) const {
				const_iterator itr =
					&myArray[( t1.x +
							   myMaxNX ) %
							 myNX][( t1.y + myMaxNY ) % myNY][( t1.z + myMaxNZ ) % myNZ];
				if( itr->second < 0 ) {
					return myEndConst;
				}
				return itr;
			}

			/// case for periodic boundary conditions
			iterator findPeriodic( const T1 &t1 ) {
				iterator itr =
					&myArray[( t1.x +
							   myMaxNX ) %
							 myNX][( t1.y + myMaxNY ) % myNY][( t1.z + myMaxNZ ) % myNZ];
				if( itr->second < 0 ) {
					return myEnd;
				}
				return itr;
			}

			T1 basisCell( const T1 &t1 ) const {
				return T1( ( t1.x + myMaxNX1 ) % myNX + myNX1,
						   ( t1.y + myMaxNY1 ) % myNY + myNY1,
						   ( t1.z + myMaxNZ1 ) % myNZ + myNZ1 );
			}

			T1 periodicCell( const T1 &t1 ) const {
				return T1( ( t1.x + myMaxNX ) % myNX,
						   ( t1.y + myMaxNY ) % myNY,
						   ( t1.z + myMaxNZ ) % myNZ );
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// New methods of class ArrayCellListStructure
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			void MaximumVolume( const float value );
			void initialize( const Vector3D &max, Vector3D cellSize );
			void updateCache();
			void uncache() {
				valid = false;
			}

			int getDimX() const {
				return myNX;
			}
			int getDimY() const {
				return myNY;
			}
			int getDimZ() const {
				return myNZ;
			}
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My data members
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			bool valid;
		private:
			TContainer myArray;
			Vector3D myCellSize;
			int myNX, myNY, myNZ;
			int myMaxNX, myMaxNY, myMaxNZ;
			int myNX1, myNY1, myNZ1;
			int myMaxNX1, myMaxNY1, myMaxNZ1;
			iterator myBegin, myEnd;
			const_iterator myBeginConst, myEndConst;
			size_type mySize;

			bool bInit;
			float mVolume, mMaxVolume;
	};
}

#endif /* ARRAYCELLSTRUCTURE_H */
