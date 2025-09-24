#ifndef PSCF_AM_ITERATOR_DARRAY_H
#define PSCF_AM_ITERATOR_DARRAY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorTmpl.h"          // base class template

namespace Pscf {

   using namespace Util;

   /**
   * Anderson mixing iterator algorithm using a DArray<double>.
   *
   * \ingroup Pscf_Iterator_Module
   */
   template <typename Iterator>
   class AmIteratorDArray 
    : public AmIteratorTmpl<Iterator, DArray<double> >
   {
   public:

      /**
      * Constructor.
      */
      AmIteratorDArray() = default;

      /**
      * Destructor.
      */
      ~AmIteratorDArray() = default;

      /**
      * Alias for base class template.
      */
      using AmTmpl = AmIteratorTmpl<Iterator, DArray<double> >;

   private:

      /**
      * Assignment for vectors of type T.
      *
      * This function must perform an assignment a = b.
      *
      * \param a  vector to be set (lhs of assignment)
      * \param b  vector value to assign (rhs of assignment)
      */
      void setEqual(DArray<double> & a, 
                    DArray<double> const & b) override;

      /**
      * Compute the inner product of two vectors.
      *
      * \param a first vector
      * \param b second vector
      */
      double dotProduct(DArray<double> const & a, 
                        DArray<double> const & b) override;


      /**
      * Return the maximum magnitude element of a vector.
      *
      * \param hist  input vector
      */
      double maxAbs(DArray<double> const & hist) override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      void subVV(DArray<double> & a, 
                 DArray<double> const & b, 
                 DArray<double> const & c) override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      void addEqVc(DArray<double> & a, 
                   DArray<double> const & b, double c) override;

   };


}
#include "AmIteratorDArray.tpp"
#endif
