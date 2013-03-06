#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/CellListEnumerator.h>
#include <protomol/topology/CubicCellManager.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Exception.h>

namespace ProtoMol {
	/**
	 * Implementation of the topology of a systems with a given boundary
	 * conditions and cell manager.
	 */

	template<class TBoundaryConditions, class TCellManager>
	class Topology : public SemiGenericTopology<TBoundaryConditions> {
		public:
			typedef CellListEnumerator<TBoundaryConditions, TCellManager> Enumerator;
			typedef TCellManager CellManager;
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Constructors, destructors, assignment
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			Topology() : SemiGenericTopology<TBoundaryConditions>() {}
			Topology( Real csf, const ExclusionType &e, const TBoundaryConditions &b,
					  const TCellManager &c, const float maxCellVolume ) : SemiGenericTopology<TBoundaryConditions>(
							  csf, e, b ), cellManager( c ) {
				cellLists.MaximumVolume( maxCellVolume );
			}

			virtual ~Topology() {};
			/// marks the the cell list as not valid any more.
			virtual void uncacheCellList() {
				cellLists.uncache();
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// New methods of class Topology
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			/// invokes an update of the cell list, if necessary
			void updateCellLists( const Vector3DBlock *positions ) const {
				if( !cellLists.valid ) {
					if( this->boundaryConditions.PERIODIC ) {
						this->min = this->boundaryConditions.getMin();
						this->max = this->boundaryConditions.getMax();
					} else {
						positions->boundingbox( this->min, this->max );
					}

					cellManager.initialize( cellLists, this->min, this->max,
											this->boundaryConditions.PERIODIC );

					Vector3D delta( this->boundaryConditions.origin() - this->min );

					/* TODO this use of CubicCellManager locks in the CellManager
					   implementation for this template.  Either the code below should be
					   changed to be more generic or the template parameter should be
					   removed.
					 */
					CubicCellManager::Cell myCell;
					CubicCellManager::CellListStructure::iterator myCellList;
					CubicCellManager::CellListStructure::iterator end = cellLists.end();

					for( int i = ( int )this->atoms.size() - 1; i >= 0; i-- ) {
						myCell =
							cellManager.findCell( delta +
												  this->boundaryConditions.minimalPosition( ( *positions )[i] ) );

						myCellList = cellLists.find( myCell );
						if( myCellList == end ) {
							// This atom is the first on its cell list, so make a new list for
							// it.
							this->atoms[i].cellListNext = -1;
							cellLists[myCell] = i;
						} else {
							this->atoms[i].cellListNext = myCellList->second;
							myCellList->second = i;
						}
					}

					cellManager.updateCache( cellLists );
				}
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From Makeable
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			virtual void getParameters( std::vector<Parameter> &parameters ) const {
				parameters.push_back( Parameter( "coulombScalingFactor",
												 Value( this->coulombScalingFactor ), 1.0 ) );
				parameters.push_back( Parameter( "exclude",
												 Value( this->exclude.getString(), ConstraintValueType::NotEmpty() ),
												 Text( std::string( "exclusion scheme (" ) +
													   ExclusionType::getPossibleValues() + std::string( ")" ) ) ) );
				this->boundaryConditions.getParameters( parameters );
				this->cellManager.getParameters( parameters );
			}

			virtual std::string getIdNoAlias() const {
				return std::string( TBoundaryConditions::keyword + TCellManager::keyword );
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From Topology
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			virtual std::string print( const Vector3DBlock *positions = NULL ) const {
				unsigned int count = 0;
				unsigned int countMult = 0;
				unsigned int countH20 = 0;
				for( unsigned int i = 0; i < this->molecules.size(); ++i ) {
					count += this->molecules[i].pairs.size();
					countH20 += ( this->molecules[i].water ? 1 : 0 );
				}

				for( unsigned int i = 0; i < this->dihedrals.size(); ++i ) {
					countMult += this->dihedrals[i].multiplicity;
				}

				std::string res;
				res +=
					"Atoms                : " + toString( this->atoms.size() ) + "\n" +
					"AtomTypes            : " + toString( this->atomTypes.size() ) + "\n" +
					"Bonds                : " + toString( this->bonds.size() ) + "\n" +
					"Angles               : " + toString( this->angles.size() ) + "\n" +
					"Dihedrals            : " + toString( this->dihedrals.size() ) + " (" +
					toString( countMult ) + ")\n" +
					"Impropers            : " + toString( this->impropers.size() ) + "\n" +
					"Molecules            : " + toString( this->molecules.size() ) + "\n" +
					"Water                : " + toString( countH20 ) + "\n" +
					"Pairs                : " + toString( count ) + "\n" +
					"Degree of freedom    : " + toString( this->degreesOfFreedom ) + "\n";
				if( positions != NULL )
					res += "Molecule pair dist's : " + toString(
							   this->checkMoleculePairDistances( *positions ) ) + "\n";
				else
					res += "Molecule pair dist's : " + toString(
							   this->minimalMolecularDistances ) + "\n";

				res +=
					"Exclusion pairs      : " +
					toString( this->exclusions.getTable().size() ) + "\n" +
					"Max exclusion dist   : " + toString( this->exclusions.getMaxDelta() ) +
					"\n" +
					"Time                 : " + toString( this->time ) + "\n" +
					"CellManager          : " + TCellManager::keyword + "\n" +
					"BoundaryConditions   : " + TBoundaryConditions::keyword +
					( this->boundaryConditions.isOrthogonal() ? " orthogonal" : "" ) + "\n";
				std::vector<Parameter> parameters;
				this->getParameters( parameters );
				for( unsigned int i = 0; i < parameters.size(); i++ )
					res +=
						getRightFill( parameters[i].keyword,
									  21 ) + ": " + parameters[i].value.getString() + "\n";

				if( positions != NULL ) {
					Real v = this->getVolume( *positions );
					if( v > Constant::EPSILON ) {
						res += "Atom density         : " +
							   toString( this->atoms.size() / v ) + "\n";
						res += "Atom/cell            : " + toString(
								   this->atoms.size() * this->cellManager.getCellVolume() / v ) + "\n";
						res += "Real cell size       : (" + toString(
								   this->cellManager.getCellSizeVector().c[0] ) + "," + toString(
								   this->cellManager.getCellSizeVector().c[1] ) + "," + toString(
								   this->cellManager.getCellSizeVector().c[2] ) + ")\n";
						res += "Cell dimension       : (" +
							   toString( this->cellLists.getDimX() ) + "," + toString(
								   this->cellLists.getDimY() ) + "," +
							   toString( this->cellLists.getDimZ() ) + ")\n";
					}
				}

				Vector3D a, b;
				this->getBoundaryConditionsBox( a, b );
				res +=
					"Simulation box       : (" + toString( a.c[0] ) + "," + toString( a.c[1] ) +
					"," + toString( a.c[2] ) + ")-" + "(" + toString( b.c[0] ) + "," +
					toString( b.c[1] ) + "," + toString( b.c[2] ) + ") " + "(" + toString( b.c[0] - a.c[0] ) +
					"," + toString( b.c[1] - a.c[1] ) + "," + toString( b.c[2] - a.c[2] ) + ")";

				if( positions != NULL ) {
					this->getBoundingbox( *positions, a, b );
					res +=
						"\nParticle             : (" + toString( a.c[0] ) + "," + toString( a.c[1] ) +
						"," + toString( a.c[2] ) + ")-" + "(" + toString( b.c[0] ) + "," +
						toString( b.c[1] ) + "," + toString( b.c[2] ) + ") " + "(" +
						toString( b.c[0] - a.c[0] ) + "," + toString( b.c[1] - a.c[1] ) + "," +
						toString( b.c[2] - a.c[2] ) + ")";

					positions->boundingbox( a, b );
					res +=
						"\nParticle extended    : (" + toString( a.c[0] ) + "," + toString( a.c[1] ) +
						"," + toString( a.c[2] ) + ")-" + "(" + toString( b.c[0] ) + "," +
						toString( b.c[1] ) + "," + toString( b.c[2] ) + ") " + "(" +
						toString( b.c[0] - a.c[0] ) + "," + toString( b.c[1] - a.c[1] ) + "," +
						toString( b.c[2] - a.c[2] ) + ")";
				}

				return res;
			}

		protected:
			virtual GenericTopology *doMake( const std::vector<Value> &values ) const {
				Real csf;
				if( !values[0].get( csf ) )
					THROW( string( " coulombScalingFactor \'" ) + values[0].getString() +
						   "\' not valid." );

				ExclusionType e( values[1].getString() );
				if( !e.valid() )
					THROW( string( " Exclusion '" ) + values[1].getString() +
						   "' not recognized, possible values are: " +
						   ExclusionType::getPossibleValues( "," ) + "." );

				// Split up parameters
				unsigned int n = TBoundaryConditions::getParameterSize();
				std::vector<Value> parmsBC( values.begin() + 2, values.begin() + n + 2 );
				std::vector<Value> parmsCM( values.begin() + n + 2, values.end() );
				
				Real maxCellVolume = 0.0f;
				if( !parmsCM[1].get( maxCellVolume ) ){
					maxCellVolume = 16.0f;
				}

				TBoundaryConditions BC = TBoundaryConditions::make( parmsBC );
				TCellManager CM = TCellManager::make( parmsCM );

				return new Topology<TBoundaryConditions, TCellManager>( csf, e, BC, CM, maxCellVolume );
			}


			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My data members
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			CellManager cellManager;
			mutable typename CellManager::CellListStructure cellLists;
	};
}
#endif /* TOPOLOGY_H */
