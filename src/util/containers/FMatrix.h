#ifndef UTIL_F_MATRIX_H
#define UTIL_F_MATRIX_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Matrix.h>

namespace Util
{

   /**
   * Fixed Size Matrix.
   *
   * The FMatrix class wraps a statically allocated 1D C array, but
   * provides access to its elements via the A(i,j) Matrix syntax.
   *
   * Template parameters M and N are the number of rows and columns
   * respectively, so that capacity1 = M and capacity2 = N.
   *
   * \ingroup Matrix_Module
   */
   template <typename Data, int M, int N>
   class FMatrix : public Matrix<Data>
   { 
   
      using Matrix<Data>::data_;
      using Matrix<Data>::capacity1_;
      using Matrix<Data>::capacity2_;
   
   public: 
   
      /**
      * Default constructor.
      */
      FMatrix();

      /**
      * Copy constructor.
      */
      FMatrix(const FMatrix<Data, M, N>& other);

      /**
      * Destructor.
      */
      ~FMatrix();

      /**
      * Assignment.
      */
      FMatrix<Data, M, N>& operator = (const FMatrix<Data, M, N>& other);

      /**
      * Serialize an FMatrix to/from an Archive.
      *
      * \param ar       archive 
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// 1D fixed size C array that holds elements of 2D array.   
      Data fixedArray_[M*N];
   
   };

   /*
   * Default constructor.
   */
   template <typename Data, int M, int N>
   FMatrix<Data,M,N>::FMatrix() 
     : Matrix<Data>()
   {
      data_ = fixedArray_;
      capacity1_ = M;
      capacity2_ = N;
   }
   
   /*
   * Copy constructor.
   */
   template <typename Data, int M, int N>
   FMatrix<Data,M,N>::FMatrix(const FMatrix<Data, M, N>& other) 
     : Matrix<Data>()
   {
      data_ = fixedArray_;
      capacity1_ = M;
      capacity2_ = N;
      for (int i = 0; i < M*N; ++i) {
         fixedArray_[i] = other.fixedArray_[i];
      }
   }
   
   /*
   * Destructor.
   */
   template <typename Data, int M, int N>
   FMatrix<Data,M,N>::~FMatrix()
   {}

   /*
   * Assignment.
   */
   template <typename Data, int M, int N>
   FMatrix<Data,M,N>& 
   FMatrix<Data,M,N>::operator = (const FMatrix<Data, M, N>& other) 
   {
      //Check for self assignment
      if (this == &other) return *this;

      // Copy elements
      for (int i = 0; i < M*N; ++i) {
         fixedArray_[i] = other.fixedArray_[i];
      }

      return *this;
   }
   
   /*
   * Serialize a FMatrix to/from an Archive.
   */
   template <class Data, int M, int N>
   template <class Archive>
   void FMatrix<Data,M,N>::serialize(Archive& ar, const unsigned int version)
   {
      for (int i = 0; i < M*N; ++i) {
         ar & fixedArray_[i];
      }
   }

} 
#endif
