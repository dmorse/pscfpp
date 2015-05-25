#ifndef UTIL_RAGGED_MATRIX_H
#define UTIL_RAGGED_MATRIX_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * A 2D array in which different rows can have different lengths.
   *
   * A RaggedMatrix object A is a two-dimensional array in which the 
   * operator A(i,j) returns a reference to element j of row i, and in
   * which different rows have different lengths. Class RaggedMatrix
   * cannot be instantiated, and functions like an abstract base class.
   *
   * The memory for a RaggedMatrix is stored in a one-dimensional C array.
   *
   * \ingroup Matrix_Module
   */
   template <typename Data>
   class RaggedMatrix
   {

   public:

      // Protected default constructor, to prohibit direct instantiation.

      // Private copy constructor, to prohibit copy construction.

      /**
      * Destructor.
      */
      virtual ~RaggedMatrix();

      /**
      * Get number of rows.
      *
      * \return Number of rows (i.e., range of first array index)
      */
      int capacity1();

      /**
      * Get number of elements in row number i.
      *
      * \param  i row index
      * \return Number of elements in row i.
      */
      int capacity2(int i);

      /**
      * Return element (i,j) of matrix by const reference.
      *
      * \param  i row index.
      * \param  j column index.
      * \return element (i, j)
      */
      const Data& operator() (int i, int j) const;

      /**
      * Return element (i,j) of matrix by reference.
      *
      * \param  i row index.
      * \param  j column index.
      * \return element (i, j)
      */
      Data& operator() (int i, int j);

   protected:

      /// One-dimensional C array of all elements.
      Data*  data_;

      /// Array of pointers to rows.
      Data**  rows_;

      /// Array containing number of elements in each row.
      int*  capacity2_;

      /// Number of rows (range of first index).
      int  capacity1_;

      /// Total number of elements.
      int  capacity_;

      /**
      * Default constructor. 
      *
      * Protected to prevent direct instantiation.
      */
      RaggedMatrix();

   private:

      /**
      * Copy constructor. 
      *
      * Private and not implemented to prohibit copying.
      */
      RaggedMatrix(const RaggedMatrix& other);

   }; 

   // Method definitions

   /**
   * Constructor (protected).
   */
   template <typename Data>
   inline RaggedMatrix<Data>::RaggedMatrix()
    : data_(0),
      rows_(0),
      capacity2_(0),
      capacity1_(0),
      capacity_(0)
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   RaggedMatrix<Data>::~RaggedMatrix()
   {}

   /*
   * Get number of rows.
   */
   template <typename Data>
   inline int RaggedMatrix<Data>::capacity1()
   { return capacity1_; }

   /*
   * Get number of columns in row i.
   */
   template <typename Data>
   inline int RaggedMatrix<Data>::capacity2(int i)
   {  return capacity2_[i]; }

   /*
   * Return element (i,j) of matrix by const reference.
   */
   template <typename Data>
   inline const Data& RaggedMatrix<Data>::operator() (int i, int j) const
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity1_);
      assert(j >= 0);
      assert(j < capacity2_[i]);
      return *(rows_[i] + j);
   }

   /*
   * Return element (i,j) of matrix by reference.
   */
   template <typename Data>
   inline Data& RaggedMatrix<Data>::operator() (int i, int j)
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity1_);
      assert(j >= 0);
      assert(j < capacity2_[i]);
      return *(rows_[i] + j);
   }

}
#endif
