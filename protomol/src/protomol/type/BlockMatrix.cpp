#include <protomol/type/BlockMatrix.h>

#include <protomol/base/Lapack.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

namespace ProtoMol
{
  // Constructors/distructors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BlockMatrix::BlockMatrix() :
      RowStart ( 0 ), Rows ( 0 ), ColumnStart ( 0 ), Columns ( 0 ),
      arraySize( 0 )
  {
  }

  BlockMatrix::BlockMatrix( const unsigned int rowStart, const unsigned int columnStart, const unsigned int rows, const unsigned int columns )
    : RowStart( rowStart ), Rows( rows ), ColumnStart( columnStart ), Columns( columns ), arraySize( rows * columns )
  {
    MyArray.resize( arraySize );
  }

  BlockMatrix::~BlockMatrix() {}

  // Methods~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockMatrix::initialize( const unsigned int rowStart, const unsigned int columnStart, const unsigned int rows, const unsigned int columns )
  {
    RowStart    = rowStart;
    ColumnStart = columnStart;

    Rows    = rows;
    Columns = columns;

    arraySize = Rows * Columns;
    MyArray.resize( arraySize );
  }

  void BlockMatrix::clear()
  {
    std::fill( MyArray.begin(), MyArray.end(), 0.0 );
  }

  // Pointer to data
  double * BlockMatrix::arrayPointer()
  {
    return &MyArray[0];
  }

  // Move block
  void BlockMatrix::blockMove( const unsigned int newRowStart, const unsigned int newColumnStart )
  {
    RowStart    = newRowStart;
    ColumnStart = newColumnStart;
  }

  // Re-size by columns
  void BlockMatrix::columnResize( const unsigned int newColumns )
  {
    Columns   = newColumns;
    arraySize = Rows * Columns;
    MyArray.resize( arraySize );
  }

  // Add 'this' with 'bm', return result
  const BlockMatrix BlockMatrix::operator+( const BlockMatrix &bm ) const
  {
    //find bounds
    const unsigned int rl = max( RowStart, bm.RowStart );
    const unsigned int rh = min( RowStart + Rows, bm.RowStart + bm.Rows );
    const unsigned int cl = max( ColumnStart, bm.ColumnStart );
    const unsigned int ch = min( ColumnStart + Columns, bm.ColumnStart + bm.Columns );

    BlockMatrix om( rl, cl, rh - rl, ch - cl );

    for ( unsigned int i = rl; i < rh; ++i ) {
      const int iTempStart  = i - om.RowStart;
      const int iOtherStart = i - bm.RowStart;
      const int iThisStart  = i -    RowStart;

      for ( unsigned int j = cl; j < ch; ++j ) {
        om.MyArray[iTempStart + ( j - om.ColumnStart ) * om.Rows] =
          MyArray[iThisStart  + ( j - ColumnStart ) * Rows] +
          bm.MyArray[iOtherStart + ( j - bm.ColumnStart ) * bm.Rows];
      }
    }

    return om;
  }

  // Add 'this' with 'bm', return result
  BlockMatrix &BlockMatrix::operator+=( const BlockMatrix &bm )
  {
    //find bounds
    const unsigned int rl = max( RowStart, bm.RowStart );
    const unsigned int rh = min( RowStart + Rows, bm.RowStart + bm.Rows );
    const unsigned int cl = max( ColumnStart, bm.ColumnStart );
    const unsigned int ch = min( ColumnStart + Columns, bm.ColumnStart + bm.Columns );

    for ( unsigned int i = rl; i < rh; ++i ) {
      const int iOtherStart = i - bm.RowStart;
      const int iThisStart  = i -    RowStart;

      for ( unsigned int j = cl; j < ch; ++j ) {
        MyArray[iThisStart + ( j - ColumnStart ) * Rows] +=
          bm.MyArray[iOtherStart + ( j - bm.ColumnStart ) * bm.Rows];
      }
    }

    return *this;
  }

  // Add 'this' with 'bm', put result in 'om'
  void BlockMatrix::add( const BlockMatrix &bm, BlockMatrix &om ) const
  {
    //find bounds
    const unsigned int rl = max( max( RowStart, bm.RowStart ), om.RowStart );
    const unsigned int rh = min( min( RowStart + Rows, bm.RowStart + bm.Rows ), om.RowStart + om.Rows );
    const unsigned int cl = max( max( ColumnStart, bm.ColumnStart ), om.ColumnStart );
    const unsigned int ch = min( min( ColumnStart + Columns, bm.ColumnStart + bm.Columns ), om.ColumnStart + om.Columns );

    for ( unsigned int i = rl; i < rh; ++i ) {
      const int iTempStart  = i - om.RowStart;
      const int iOtherStart = i - bm.RowStart;
      const int iThisStart  = i -    RowStart;

      for ( unsigned int j = cl; j < ch; ++j ) {
        om.MyArray[iTempStart + ( j - om.ColumnStart ) * om.Rows] =
          MyArray[iThisStart  + ( j - ColumnStart ) * Rows] +
          bm.MyArray[iOtherStart + ( j - bm.ColumnStart ) * bm.Rows];
      }
    }
  }

  // Multiply 'this' with 'bm', return result: TEST A, C
  const BlockMatrix BlockMatrix::operator*( const BlockMatrix &bm ) const
  {
    const unsigned int kl = max( ColumnStart, bm.RowStart );
    const unsigned int kh = min( ColumnStart + Columns, bm.RowStart + bm.Rows );

    BlockMatrix om( RowStart, bm.ColumnStart, Rows, bm.Columns );

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"N"; char *transB = ( char * )"N";
		int m = Rows; int n = bm.Columns; int k = kh - kl;
		int lda = Rows; int ldb = bm.Rows; int ldc = om.Rows;
		double alpha = 1.0; double beta = 0.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-ColumnStart )*Rows], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
             &ldb, &beta, &om.MyArray[( RowStart - om.RowStart ) + ( bm.ColumnStart - om.ColumnStart ) * om.Rows], &ldc );
	}else{
		const unsigned int rowLength = RowStart + Rows;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = RowStart; i < rowLength; ++i ) {
		  const int iThisRowStart  = i - RowStart;
		  const int iOtherRowStart = i - om.RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;

			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[ iThisRowStart + ( k - ColumnStart ) * Rows] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}

			om.MyArray[ iOtherRowStart + ( j - om.ColumnStart ) * om.Rows] = tempf; //.at( changed to [] for efficiency
		  }
		}
	}

    return om;
  }

  // Multiply 'this' with 'bm', put result in 'om': TEST A, B, C
  void BlockMatrix::product( const BlockMatrix &bm, BlockMatrix &om ) const
  {
    const unsigned int kl = max( ColumnStart, bm.RowStart );
    const unsigned int kh = min( ColumnStart + Columns, bm.RowStart + bm.Rows );

    if ( om.RowStart > RowStart || om.ColumnStart > bm.ColumnStart ||
         om.RowStart + om.Rows < RowStart + Rows || om.ColumnStart + om.Columns < bm.ColumnStart + bm.Columns )
      Report::report << Report::error << "[BlockMatrix::product] Target Block Matrix too small." << Report::endr;

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"N"; char *transB = ( char * )"N";
		int m = Rows; int n = bm.Columns; int k = kh - kl;
		int lda = Rows; int ldb = bm.Rows; int ldc = om.Rows;
		double alpha = 1.0; double beta = 0.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-ColumnStart )*Rows], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
             &ldb, &beta, &om.MyArray[( RowStart - om.RowStart ) + ( bm.ColumnStart - om.ColumnStart ) * om.Rows], &ldc );
	}else{
		const unsigned int rowLength = RowStart + Rows;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = RowStart; i < rowLength; ++i ) {
		  const int iThisRowStart  = i - RowStart;
		  const int iOtherRowStart = i - om.RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;

			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[ iThisRowStart + ( k - ColumnStart ) * Rows] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}

			om.MyArray[ iOtherRowStart + ( j - om.ColumnStart ) * om.Rows] = tempf; //.at( changed to [] for efficiency
		  }
		}
	}
  }

  // Multiply 'this' with 'bm', put result in double array 'om': TEST A, B, C
  void BlockMatrix::productToArray( const BlockMatrix &bm,
                                    double *om_MyArray, unsigned int om_RowStart, unsigned int om_ColumnStart,
                                    unsigned int om_Rows, unsigned int om_Columns ) const
  {
    const unsigned int kl = max( ColumnStart, bm.RowStart );
    const unsigned int kh = min( ColumnStart + Columns, bm.RowStart + bm.Rows );

    if ( om_RowStart > RowStart || om_ColumnStart > bm.ColumnStart ||
         om_RowStart + om_Rows < RowStart + Rows || om_ColumnStart + om_Columns < bm.ColumnStart + bm.Columns )
      Report::report << Report::error << "[BlockMatrix::product] Target Block Matrix too small." << Report::endr;

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"N"; char *transB = ( char * )"N";
		int m = Rows; int n = bm.Columns; int k = kh - kl;
		int lda = Rows; int ldb = bm.Rows; int ldc = om_Rows;
		double alpha = 1.0; double beta = 0.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-ColumnStart )*Rows], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
             &ldb, &beta, &om_MyArray[( RowStart - om_RowStart ) + ( bm.ColumnStart - om_ColumnStart ) * om_Rows], &ldc );
	}else{
		const unsigned int rowLength = RowStart + Rows;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = RowStart; i < rowLength; ++i ) {
		  const int iThisRowStart  = i - RowStart;
		  const int iOtherRowStart = i - om_RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;

			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[ iThisRowStart + ( k - ColumnStart ) * Rows] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}

			om_MyArray[ iOtherRowStart + ( j - om_ColumnStart ) * om_Rows] = tempf; //.at( changed to [] for efficiency
		  }
		}
	}
  }

  // Multiply 'this' with 'bm', sum result in 'om': TEST A, B, C
  void BlockMatrix::sumProduct( const BlockMatrix &bm, BlockMatrix &om ) const
  {
    const unsigned int kl = max( ColumnStart, bm.RowStart );
    const unsigned int kh = min( ColumnStart + Columns, bm.RowStart + bm.Rows );

    if ( om.RowStart > RowStart || om.ColumnStart > bm.ColumnStart ||
         om.RowStart + om.Rows < RowStart + Rows || om.ColumnStart + om.Columns < bm.ColumnStart + bm.Columns )
      Report::report << Report::error << "[BlockMatrix::product] Target Block Matrix too small." << Report::endr;

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"N"; char *transB = ( char * )"N";
		int m = Rows; int n = bm.Columns; int k = kh - kl;
		int lda = Rows; int ldb = bm.Rows; int ldc = om.Rows;
		double alpha = 1.0; double beta = 1.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-ColumnStart )*Rows], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
             &ldb, &beta, &om.MyArray[( RowStart - om.RowStart ) + ( bm.ColumnStart - om.ColumnStart ) * om.Rows], &ldc );
	}else{
		const unsigned int rowLength = RowStart + Rows;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = RowStart; i < rowLength; ++i ) {
		  const int iThisRowStart  = i - RowStart;
		  const int iOtherRowStart = i - om.RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;

			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[ iThisRowStart + ( k - ColumnStart ) * Rows] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}

			om.MyArray[ iOtherRowStart + ( j - om.ColumnStart ) * om.Rows] += tempf; //.at( changed to [] for efficiency
		  }
		}
	}
  }

  // Multiply transpose of 'this' with 'bm', put result in 'om': TEST A, B, C
  void BlockMatrix::transposeProduct( const BlockMatrix &bm, BlockMatrix &om ) const
  {
    const unsigned int kl = max( RowStart, bm.RowStart );
    const unsigned int kh = min( RowStart + Rows, bm.RowStart + bm.Rows );

    if ( om.RowStart > ColumnStart || om.ColumnStart > bm.ColumnStart ||
         om.RowStart + om.Rows < ColumnStart + Columns || om.ColumnStart + om.Columns < bm.ColumnStart + bm.Columns )
      Report::report << Report::error << "[BlockMatrix::transposeProduct] Target Block Matrix too small." << Report::endr;

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"T"; char *transB = ( char * )"N";
		int m = Columns; int n = bm.Columns; int k = kh - kl;//
		int lda = Rows; int ldb = bm.Rows; int ldc = om.Rows;//
		double alpha = 1.0; double beta = 0.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-RowStart )], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
				 &ldb, &beta, &om.MyArray[( ColumnStart - om.RowStart ) + ( bm.ColumnStart - om.ColumnStart ) * om.Rows], &ldc );
	}else{
		const unsigned int columnLength = ColumnStart + Columns;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = ColumnStart; i < columnLength; ++i ) {
		  const int iColumnStart = i - ColumnStart;
		  const int iRowStart    = i - om.RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;
			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[ iColumnStart * Rows + ( k - RowStart ) ] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}
			om.MyArray.at( iRowStart + ( j - om.ColumnStart ) * om.Rows ) = tempf; //.at( changed to [] for efficiency
		  }
		}
	}
  }

  // Multiply transpose of 'this' with 'bm', return result: TEST C
  const BlockMatrix BlockMatrix::operator/( const BlockMatrix &bm ) const
  {
    const unsigned int kl = max( RowStart, bm.RowStart );
    const unsigned int kh = min( RowStart + Rows, bm.RowStart + bm.Rows );

    BlockMatrix om( ColumnStart, bm.ColumnStart, Columns, bm.Columns );

	if( Lapack::isEnabled() ){
		char *transA = ( char * )"T"; char *transB = ( char * )"N";
		int m = Columns; int n = bm.Columns; int k = kh - kl;
		int lda = Rows; int ldb = bm.Rows; int ldc = om.Rows;
		double alpha = 1.0; double beta = 0.0;

		Lapack::dgemm( transA, transB, &m, &n, &k, &alpha, ( double * )&MyArray[( kl-RowStart )], &lda, ( double * )&bm.MyArray[kl - bm.RowStart],
				 &ldb, &beta, &om.MyArray[( RowStart - om.RowStart ) + ( bm.ColumnStart - om.ColumnStart ) * om.Rows], &ldc );
	}else{
		const unsigned int columnLength = ColumnStart + Columns;
		const unsigned int otherColumnLength = bm.ColumnStart + bm.Columns;

		for ( unsigned int i = ColumnStart; i < columnLength; ++i ) {
		  const int iColumnStart = i - ColumnStart;
		  const int iRowStart    = i - om.RowStart;

		  for ( unsigned int j = bm.ColumnStart; j < otherColumnLength; ++j ) {
			Real tempf = 0.0;

			for ( unsigned int k = kl; k < kh; ++k ) {
			  tempf += MyArray[iColumnStart * Rows + ( k - RowStart ) ] *
					   bm.MyArray[( k - bm.RowStart ) + ( j - bm.ColumnStart ) * bm.Rows];
			}

			om.MyArray[iRowStart + ( j - om.ColumnStart ) * om.Rows] = tempf;
		  }
		}
	}

    return om;
  }

  // Get sub-matrix: TEST A
  const BlockMatrix BlockMatrix::subMatrix( const unsigned int atRow, const unsigned int atColumn,
      const unsigned int getRows, const unsigned int getColumns ) const
  {
    //create matrix
    BlockMatrix om( atRow, atColumn, getRows, getColumns );
    om.clear();

    const unsigned int atRowSize = atRow + getRows;
    const unsigned int atColSize = atColumn + getColumns;

    if ( atRow >= RowStart && atRowSize <= RowStart + Rows &&
         atColumn >= ColumnStart && atColSize <= ColumnStart + Columns ) {

      for ( unsigned int i = atRow; i < atRowSize; ++i ) {
        const int iTempStart  = i - om.RowStart;
        const int iThisStart  = i -    RowStart;

        for ( unsigned int j = atColumn; j < atColSize; ++j ) {
          om.MyArray.at( iTempStart + ( j - om.ColumnStart ) * om.Rows ) =
            MyArray[iThisStart + ( j - ColumnStart ) * Rows];
        }
      }
    }
    return om;
  }

  // Operators~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Double Index access
  double & BlockMatrix::operator()(const unsigned int rowIndex,
                                   const unsigned int colIndex )
  {
    if ( rowIndex < RowStart || rowIndex > RowStart + Rows || colIndex < ColumnStart || colIndex > ColumnStart + Columns ) {
      Report::report << Report::error << "[BlockMatrix::operator(,)] Index out of range." << Report::endr;
    }

    return MyArray[( rowIndex - RowStart ) + ( colIndex - ColumnStart ) * Rows];
  }

  // Index access
  double & BlockMatrix::operator[](const unsigned int index )
  {
    if ( index > arraySize ) {
      Report::report << Report::error << "[BlockMatrix::operator[]] Index out of range." << Report::endr;
    }

    return MyArray[index];
  }

  // Equals
  BlockMatrix & BlockMatrix::operator=( const BlockMatrix & bm )
  {
    arraySize = bm.arraySize;

    MyArray.resize( arraySize );

    std::copy ( bm.MyArray.begin(), bm.MyArray.end(), MyArray.begin() );

    RowStart    = bm.RowStart;
    ColumnStart = bm.ColumnStart;

    Rows    = bm.Rows;
    Columns = bm.Columns;

    return *this;
  }
}
