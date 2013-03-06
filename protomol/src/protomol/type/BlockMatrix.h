/*  -*- c++ -*-  */
#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <vector>
#include <algorithm>
#include <protomol/base/Report.h>
#include <protomol/type/Real.h>

using namespace std;

namespace ProtoMol
{
  //_________________________________________________________________ BlockMatrix
  /**
   * Block Matrix object
   */

  class BlockMatrix
  {

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      unsigned int RowStart, Rows, ColumnStart, Columns;
      vector<double> MyArray;
      unsigned int arraySize;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      BlockMatrix();

      BlockMatrix( const unsigned int rowStart, const unsigned int columnStart, const unsigned int rows, const unsigned int columns );

      ~BlockMatrix();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class Hessian
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      void initialize( const unsigned int rowStart, const unsigned int columnStart, const unsigned int rows, const unsigned int columns );
      void clear();

      // Re-size by columns
      void columnResize( const unsigned int newColumns );

      // Move block
      void blockMove( const unsigned int newRowStart, const unsigned int newColumnStart );

      // Double Index access
      double &operator()( const unsigned int rowIndex, const unsigned int colIndex );

      // Index access
      double &operator[]( const unsigned int index );

      // Equals
      BlockMatrix &operator=( const BlockMatrix &bm );

      // Add 'this' with 'bm', return result
      const BlockMatrix operator+( const BlockMatrix &bm ) const;

      // Add 'this' with 'bm', put result in this
      BlockMatrix & operator+=( const BlockMatrix &bm );

      // Add 'this' with 'bm', put result in 'om'
      void add( const BlockMatrix &bm, BlockMatrix &om ) const;

      // Multiply 'this' with 'bm', return result
      const BlockMatrix operator*( const BlockMatrix &bm ) const;

      // Multiply 'this' with 'bm', put result in 'om'
      void product( const BlockMatrix &bm, BlockMatrix &om ) const;

      // Multiply 'this' with 'bm', put result in double array 'om': TEST A, B, C
      void productToArray( const BlockMatrix &bm,
                           double *om_MyArray, unsigned int om_RowStart, unsigned int om_ColumnStart,
                           unsigned int om_Rows, unsigned int om_Columns ) const;

      // Multiply 'this' with 'bm', sum result in 'om'
      void sumProduct( const BlockMatrix &bm, BlockMatrix &om ) const;

      // Multiply transpose of 'this' with 'bm', put result in 'om'
      void transposeProduct( const BlockMatrix &bm, BlockMatrix &om ) const;

      // Multiply transpose of 'this' with 'bm', return result: TEST
      const BlockMatrix operator/( const BlockMatrix &bm ) const;

      // Get sub-matrix
      const BlockMatrix subMatrix( const unsigned int atRow, const unsigned int atColumn, const unsigned int getRows, const unsigned int getColumns ) const;

      // Pointer to data
      double * arrayPointer();

  };

}

#endif

