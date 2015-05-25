#ifndef UTIL_TENSOR_H
#define UTIL_TENSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/space/Dimension.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif

#include <iostream>
#include <cmath>

namespace Util
{

   class Vector;

   /**
   * A Tensor represents a Cartesian tensor.
   *
   * \ingroup Space_Module
   */
   class Tensor
   {

   public:

      /// \name Constructors
      //@{ 
      
      /**
      * Default constructor
      */
      Tensor();

      /**
      * Copy constructor
      */
      Tensor(const Tensor& t);

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param scalar initial value for all elements.
      */
      explicit Tensor(double scalar);

      /**
      * Construct Tensor from double [][Dimension] 2D C array.
      *
      * \param a 2D array a[Dimension][Dimension]
      */
      explicit Tensor(const double a[][Dimension]);

      //@}
      
      /// \name Assignment
      //@{
      
      /**
      * Copy assignment.
      *
      * \param t Tensor to assign.
      */
      Tensor& operator=(const Tensor& t);

      /**
      * Assignment from C double [Dimension][Dimension] 2D array.
      *
      * \param a 2D array a[Dimension][Dimension]
      */
      Tensor& operator=(const double a[][Dimension]);

      //@}
      /// \name Arithmetic Assignment
      //@{
      
      /**
      * Add tensor dt to this tensor.
      *
      * Upon return, *this = *this + dt.
      *
      * \param dt  tensor increment (input)
      */
      void operator+=(const Tensor& dt);

      /**
      * Subtract tensor dt from this tensor.
      *
      * Upon return, *this = *this - dt.
      *
      * \param dt   tensor increment (input)
      */
      void operator-=(const Tensor& dt);

      /**
      * Multiply this tensor by scalar s.
      *
      * Upon return, *this = (*this)*s.
      *
      * \param s  scalar multiplier
      */
      void operator*=(double s);

      /**
      * Divide this tensor by scalar s.
      *
      * Upon return, *this = (*this)/s.
      *
      * \param s  scalar divisor (input)
      */
      void operator/=(double s);

      //@}
      /// \name Array Subscript 
      //@{
      
      /**
      * Return one element by value.
      *
      * \param i  row element index
      * \param j  column element index
      * \return element (i, j) of the tensor
      */
      const double& operator ()(int i, int j) const;

      /**
      * Return one element by non-const reference.
      *
      * \param i  row element index
      * \param j  column element index
      * \return  element i of the tensor
      */
      double& operator () (int i, int j);

      //@}
      /// \name Tensor valued functions (result assigned to invoking object)
      //@{
      
      /**
      * Set all elements of this tensor to zero.
      *
      * \return reference to this tensor
      */
      Tensor& zero();

      /**
      * Set row i of this Tensor to elements of Vector r.
      *
      * \return reference to this tensor
      */
      Tensor& setRow(int i, const Vector& r);

      /**
      * Set column i of this Tensor to elements of Vector r.
      *
      * \return reference to this tensor
      */
      Tensor& setColumn(int i, const Vector& r);

      /**
      * Add tensors t1 and t2.
      *
      * Upon return, *this = t1 + t2.
      *
      * \param t1  tensor
      * \param t2  tensor
      * \return reference to this tensor
      */
      Tensor& add(const Tensor& t1, const Tensor& t2);

      /**
      * Subtract tensor t2 from t1.
      *
      * Upon return, *this == t1 - t2.
      *
      * \param t1  tensor (input)
      * \param t2  tensor (input)
      * \return reference to this tensor
      */
      Tensor& subtract(const Tensor& t1, const Tensor& t2);

      /**
      * Multiply a tensor t by a scalar s.
      *
      * Upon return, *this == v*s.
      *
      * \param t  tensor factor
      * \param s  scalar factor
      * \return reference to this tensor
      */
      Tensor& multiply(const Tensor& t, double s);

      /**
      * Divide a Tensor t by a scalar s.
      *
      * Upon return, *this = v/s;
      *
      * \param t  tensor input
      * \param s  scalar denominator
      * \return reference to this tensor
      */
      Tensor& divide(const Tensor& t, double s);

      /**
      * Compute transpose of a tensor.
      *
      * Upon return, *this is the transpose of t
      *
      * \param  t  input tensor
      * \return reference to this tensor
      */
      Tensor& transpose(const Tensor& t);

      /**
      * Transpose this tensor.
      *
      * Upon return, *this is transposed.
      *
      * \return reference to this tensor
      */
      Tensor& transpose();

      /**
      * Compute symmetric part of a tensor t.
      *
      * Upon return, *this = [t + t.transpose()]/2
      *
      * \param  t  tensor input
      * \return reference to this tensor
      */
      Tensor& symmetrize(const Tensor& t);

      /**
      * Symmetrize this tensor.
      *
      * Upon return, this is symmetrized, equal to half the
      * sum of the original tensor and its transpose.
      *
      * \return reference to this tensor
      */
      Tensor& symmetrize();

      /**
      * Create dyad of two vectors.
      *
      * Upon return, *this equals the dyad v1 ^ v2.
      * Equivalently: (*this)(i , j) == v1[i]*v2[j]
      *
      * \param  v1 vector input
      * \param  v2 vector input
      * \return reference to this tensor
      */
      Tensor& dyad(const Vector& v1, const Vector& v2);

      //@}
      /// \name Miscellaneous 
      //@{ 

      /**
      * Return the trace of this tensor.
      */
      double trace() const; 

      /**
      * Serialize this to/from an archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      //@}
      /// \name Static Members
      //@{ 
     
      /** 
      * Constant Tensor with all zero elements.
      */
      static const Tensor Zero;

      /**
      * Call to guarantee initialization of Tensor::Zero tensor.
      */
      static void initStatic();

      #ifdef UTIL_MPI
      /**
      * Commit MPI datatype MpiTraits<Tensor>::type.
      */
      static void commitMpiType();
      #endif

      //@}

   private:

      /// Width of field in stream IO
      static const int Width = 13;

      /// Precision in stream IO of Tensor coordinates
      static const int Precision = 5;

      /// Elements of the tensor.
      double elem_[DimensionSq];

   //friends:

      friend bool operator == (const Tensor& t1, const Tensor& t2);

      friend bool operator == (const Tensor& t1, const double t2[][Dimension]);

      friend std::istream& operator >> (std::istream& in, Tensor &tensor);

      friend std::ostream& operator << (std::ostream& out, const Tensor &tensor);

   };

   /// Equality for Tensors.
   bool operator == (const Tensor& t1, const Tensor& t2);

   /// Equality of Tensor and 2D C array.
   bool operator == (const Tensor& t1, const double t2[][Dimension]);

   /// Equality of C array and Tensor.
   bool operator == (const double t1[][Dimension], const Tensor& t2);

   /// Inequality of two Tensors.
   bool operator != (const Tensor& t1, const Tensor& t2);

   /// Inequality of Tensor and C array.
   bool operator != (const Tensor& t1, const double t2[][Dimension]);

   /// Inequality of C array and Tensor.
   bool operator != (const double t1[][Dimension], const Tensor& t2);

   /**
   * istream extractor for a Tensor.
   *
   * Input elements of a tensor from stream, without line breaks.
   *
   * \param in      input stream
   * \param tensor  Tensor to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Tensor &tensor);

   /**
   * ostream inserter for a Tensor.
   *
   * Output elements of a tensor to stream, without line breaks.
   * \param  out     output stream
   * \param  tensor  Tensor to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Tensor &tensor);

   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<Tensor>.
   */
   template <>
   class MpiTraits<Tensor>
   {  
   public:  
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };
   #endif

   // Inline methods

   /*
   * Default constructor
   */
   inline
   Tensor::Tensor()
   {}

   /*
   * Copy constructor
   */
   inline
   Tensor::Tensor(const Tensor& t)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t.elem_[i];
      }
   }

   /*
   * Constructor, initialize all elements to a scalar value.
   */
   inline
   Tensor::Tensor(double scalar)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = scalar;
      }
   }

   /*
   * Construct Tensor from C double [][Dimension] array.
   */
   inline
   Tensor::Tensor(const double a[][Dimension])
   {
      for (int i = 0; i < Dimension; ++i) {
         for (int j = 0; j < Dimension; ++j) {
            elem_[i*Dimension + j] = a[i][j];
         }
      }
   }

   /*
   * Set all elements of this tensor to zero.
   */
   inline
   Tensor& Tensor::zero()
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = 0.0;
      }
      return *this;
   }

   /*
   * Assignment.
   */
   inline
   Tensor& Tensor::operator=(const Tensor& t)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t.elem_[i];
      }
      return *this;
   }

   /*
   * Assignment from C double [][Dimension] array.
   */
   inline
   Tensor& Tensor::operator=(const double a[][Dimension])
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            elem_[i*Dimension + j] = a[i][j];
         }
      }
      return *this;
   }

   /*
   * Add tensor dt to this tensor.
   */
   inline
   void Tensor::operator+=(const Tensor& dt)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] += dt.elem_[i];
      }
   }

   /*
   * Subtract tensor dt from this tensor.
   */
   inline
   void Tensor::operator-=(const Tensor& dt)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] -= dt.elem_[i]; 
      }
   }

   /*
   * Multiply this tensor by scalar s.
   */
   inline
   void Tensor::operator*=(double s)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] *= s; 
      }
   }

   /*
   * Divide this tensor by scalar s.
   */
   inline
   void Tensor::operator/=(double s)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] /= s; 
      }
   }

   /*
   * Return one element by value.
   */
   inline
   const double& Tensor::operator()(int i, int j) const
   {
      assert(i >=0);  
      assert(i < Dimension);  
      assert(j >=0);  
      assert(j < Dimension);  
      return elem_[i*Dimension + j]; 
   }

   /*
   * Return a reference to one element of the tensor.
   */
   inline
   double& Tensor::operator()(int i, int j)
   {  
      assert(i >=0);  
      assert(i < Dimension);  
      assert(j >=0);  
      assert(j < Dimension);  
      return elem_[i*Dimension + j]; 
   }

   /*
   * Add tensors t1 and t2.
   *
   * Upon return, *this = t1 + t2.
   */
   inline
   Tensor& Tensor::add(const Tensor& t1, const Tensor& t2)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t1.elem_[i] + t2.elem_[i];
      }
      return *this;
   }

   /*
   * Subtract tensor t2 from t1.
   *
   * Upon return, *this = t1 - t2.
   */
   inline
   Tensor& Tensor::subtract(const Tensor& t1, const Tensor& t2)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t1.elem_[i] - t2.elem_[i];
      }
      return *this;
   }

   /*
   * Multiply a tensor t by a scalar s.
   *
   * Upon return, *this = t*s.
   */
   inline
   Tensor& Tensor::multiply(const Tensor& t, double s)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t.elem_[i]*s;
      }
      return *this;
   }

   /*
   * Divide tensor t by scalar s.
   *
   * Upon return, *this = t/s;
   */
   inline
   Tensor& Tensor::divide(const Tensor& t, double s)
   {
      for (int i = 0; i < DimensionSq; ++i) {
         elem_[i] = t.elem_[i]/s;
      }
      return *this;
   }

   /*
   * Return trace of a tensor.
   */
   inline
   double Tensor::trace() const
   {
      double trace = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         trace += (*this)(i, i);
      }
      return trace;
   }

   /*
   * Compute the transpose of a tensor.
   *
   * Upon return, this tensor is the transpose of t.
   */
   inline
   Tensor& Tensor::transpose(const Tensor& t)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         (*this)(i, i) = t(i, i);
      }
      for (i = 1; i < Dimension; ++i) {
         for (j = 0; j < i; ++j) {
            (*this)(i, j) = t(j, i);  
            (*this)(j, i) = t(i, j);
         }
      }
      return *this;
   }

   /*
   * Transpose this tensor.
   *
   * Upon return, this tensor is transposed.
   */
   inline
   Tensor& Tensor::transpose()
   {
      double save;
      int i, j;
      for (i = 1; i < Dimension; ++i) {
         for (j = 0; j < i; ++j) {
            save = (*this)(i, j);
            (*this)(i, j) = (*this)(j, i);  
            (*this)(j, i) = save;
         }
      }
      return *this;
   }

   /*
   * Compute symmetric part of a tensor.
   *
   * Upon return, this is symmetric part of t:
   * *this = 0.5*(t + t.transpose());
   */
   inline
   Tensor& Tensor::symmetrize(const Tensor& t)
   {
      double ave;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         (*this)(i, i) = t(i, i);
      }
      for (i = 1; i < Dimension; ++i) {
         for (j = 0; j < i; ++j) {
            ave = 0.5*( t(i, j) + t(j, i) );
            (*this)(i, j) = ave;  
            (*this)(j, i) = ave;
         }
      }
      return *this;
   }

   /*
   * Symmetrize this tensor.
   *
   * Upon return, t is symmetrized:
   */
   inline
   Tensor& Tensor::symmetrize()
   {
      double ave;
      int i, j;
      for (i = 1; i < Dimension; ++i) {
         for (j = 0; j < i; ++j) {
            ave = 0.5*( (*this)(i, j) + (*this)(j, i) );
            (*this)(i, j) = ave;  
            (*this)(j, i) = ave;
         }
      }
      return *this;
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   inline 
   void Tensor::serialize(Archive& ar, const unsigned int version)
   { 
      for (int i = 0; i < DimensionSq; ++i) {
         ar & elem_[i];
      }
   }

}
 
#include <util/space/Vector.h>
namespace Util
{

   /*
   * Set row i of the tensor to elements of vector r.
   */
   inline
   Tensor& Tensor::setRow(int i, const Vector& r)
   {
      for (int j = 0; j < Dimension; j++) {
         elem_[i*Dimension+j] = r[j];
      }
      return *this;
   }

   /*
   * Set column i of the tensor to elements of vector r.
   */
   inline
   Tensor& Tensor::setColumn(int j, const Vector& r)
   {
      for (int i = 0; i < Dimension; i++) {
         elem_[i*Dimension+j] = r[i];
      }
      return *this;
   }

   /*
   * Multiply two vectors to create a dyadic tensor.
   *
   * Upon return, *this = t*s.
   */
   inline
   Tensor& Tensor::dyad(const Vector& v1, const Vector& v2)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            elem_[i*Dimension + j] = v1[i]*v2[j];
         }
      }
      return *this;
   }

} 
#endif
