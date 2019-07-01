#ifndef PSCF_UNIT_CELL_BASE_H
#define PSCF_UNIT_CELL_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/FArray.h>
#include <util/containers/FMatrix.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>
#include <util/math/Constants.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Base class template for a crystallographic unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   class UnitCellBase 
   {
   public:

      /**
      * Constructor.
      */
      UnitCellBase();
   
      /**
      * Destructor.
      */
      ~UnitCellBase()
      {}
   
      /** 
      * Set the parameters of unit cell.
      */  
      void SetParams(double val, int m);

      /**
      * Compute square magnitude of reciprocal basis vector.
      */
      virtual double ksq(IntVec<D> const & k) const;
   
      /**
      * Get the number of Parameters in the unit cell.
      */
      const int nParams() const;

      /** 
      * Get the parameters of unit cell.
      */  
      FArray<double, 6> params();

      /**
      * Get Bravais basis vector i, denoted by a_i.
      */
      const RealVec<D>& rBasis(int i) const;
   
      /**
      * Get reciprocal basis vector i, denoted by b_i.
      */
      const RealVec<D>& kBasis(int i) const;

      /**
      * Get component j of derivative of rBasis vector ai w/respect to k.
      */
      const double drBasis(int k, int i, int j) const;

      /**
      * Get component j of derivative of kBasis vector bi w/respect to k.
      */
      const double dkBasis(int k, int i, int j) const;
 
      /** 
      * Get the derivative of dot product ri.rj with respect to parameter k.
      */  
      const double drrBasis(int k, int i, int j) const;
  
      /**
      * Get the derivative of dot product bi.bj with respect to parameter k.
      */
      const double dkkBasis(int k, int i, int j) const;

   protected:
 
      /**
      * Array of Bravais lattice basis vectors.
      */ 
      FArray<RealVec<D>, D> rBasis_;

      /**
      * Array of reciprocal lattice basis vectors.
      */ 
      FArray<RealVec<D>, D> kBasis_;

      /**
      * Array of derivatives of rBasis.
      *
      * Element drBasis_[k](i,j) is derivative with respect to
      * parameter k of component j of Bravais basis vector i.
      */ 
      FArray<FMatrix<double, D, D>, 6> drBasis_;
      
      /**
      * Array of derivatives of kBasis
      *
      * Element dkBasis_[k](i,j) is derivative with respect to
      * parameter k of component j of reciprocal basis vector i.
      */
      FArray<FMatrix<double, D, D>, 6> dkBasis_;

      /**
      * Array of derivatives of a_i.a_j
      *
      * Element drrBasis_[k](i,j) is derivative with respect 
      * to parameter k of the dot product (a_i.a_j) of Bravais 
      * lattice basis vectors a_i and a_j.
      */
      FArray<FMatrix<double, D, D>, 6> drrBasis_;

      /**
      * Array of derivatives of b_i.b_j
      * 
      * Element dkkBasis_[k](i,j) is derivative with respect 
      * to parameter k of the dot product (b_i.b_j) of reciprocal
      * lattice basis vectors b_i and b_j.
      */
      FArray<FMatrix<double, D, D>, 6> dkkBasis_;

      /**
      * Parameters used to describe the unit cell.
      */      
      FArray<double, 6> parameters_;
   
      /**
      * Number of parameters required to specify unit cell.
      */
      int nParameter_;

      // Protected functions (used in implementation).

      /**
      * Compute all private data, given latticeSystem and parameters.
      *
      * Calls initializeToZero, setBasis, computeDerivatives internally.
      */
      void setLattice();

   private:

      /**
      * Initialize all arrays to zero.
      *
      * Sets all elements of the following arrays to zero:
      * rBasis_, kBasis_, drBasis_, dkBasis, drrBasis_ and dkkBasis_.
      */
      void initializeToZero();

      /**
      * Set values of rBasis, kBasis, and drBasis.
      *
      * Invoke initializeToZero before this, computeDerivatives after.
      *
      * \pre Lattice system, nParam and parameters must be set.
      */
      virtual void setBasis() = 0;

      /**
      * Compute values of dkBasis_, drrBasis_, and dkkBasis_.
      */
      void computeDerivatives();

   };

   /*
   * Get the number of Parameters in the unit cell.
   */
   template <int D>
   inline
   const int UnitCellBase<D>::nParams() const
   {  return nParameter_;  }

   /*
   * Get the unit cell parameters.
   */
   template <int D>
   inline
   FArray<double, 6> UnitCellBase<D>::params()
   {  return parameters_;  }

   /*
   * Get Bravais basis vector i.
   */
   template <int D>
   const RealVec<D>& UnitCellBase<D>::rBasis(int i) const
   {  return rBasis_[i];  }

   /*
   * Get reciprocal basis vector i.
   */
   template <int D>
   inline
   const RealVec<D>& UnitCellBase<D>::kBasis(int i) const
   {  return kBasis_[i];  }

   /*
   * Get component j of derivative of r basis vector i w/respect to param k.
   */
   template <int D>
   inline
   const double UnitCellBase<D>::drBasis(int k, int i, int j) const
   {  return drBasis_[k](i,j);  }

   /*
   * Get component j of derivative of r basis vector i w/respect to param k.
   */
   template <int D>
   inline
   const double UnitCellBase<D>::dkBasis(int k, int i, int j) const
   {  return dkBasis_[k](i, j);  }

   /*
   * Get the derivative of ki*kj with respect to parameter k.
   */
   template <int D>
   inline
   const double UnitCellBase<D>::dkkBasis(int k, int i, int j) const
   {  return dkkBasis_[k](i, j);  }

   /*
   * Get the derivative of ri*rj with respect to parameter k.
   */
   template <int D>
   inline
   const double UnitCellBase<D>::drrBasis(int k, int i, int j) const
   {  return drrBasis_[k](i, j);  }

  /*
   * Set the parameters in the unit cell.
   */
   template <int D>
   inline
   void UnitCellBase<D>::SetParams(double val, int m)
   {  parameters_ [m] = val;}

   /*
   * Constructor.
   */
   template <int D>
   UnitCellBase<D>::UnitCellBase()
    : nParameter_(0)
   {}

   /*
   * Get square magnitude of reciprocal basis vector.
   */
   template <int D>
   double UnitCellBase<D>::ksq(IntVec<D> const & k) const
   {
      RealVec<D> g(0.0);
      RealVec<D> p;
      for (int i = 0; i < D; ++i) {
         p.multiply(kBasis_[i], k[i]);
         g += p;
      }
      double value = 0.0;
      for (int i = 0; i < D; ++i) {
         value += g[i]*g[i];
      }
      return value;
   }

   // Protected member functions

   /*
   * Initialize internal arrays to zero.
   */
   template <int D>
   void UnitCellBase<D>::initializeToZero()
   {
      // Initialize all elements to zero
      int i, j, k;
      for (i = 0; i < D; ++i) { 
         for (j = 0; j < D; ++j) { 
            rBasis_[i][j] = 0.0;
            kBasis_[i][j] = 0.0;
         }
      }
      for (k = 0; k < 6; ++k){
         for (i = 0; i < D; ++i) { 
            for (j = 0; j < D; ++j) { 
               drBasis_[k](i,j) = 0.0;
               dkBasis_[k](i,j) = 0.0;
               drrBasis_[k](i,j) = 0.0;
               dkkBasis_[k](i,j) = 0.0;
            }
         }
      }
   }

   /*
   * Compute quantities involving derivatives.
   */
   template <int D>
   void UnitCellBase<D>::computeDerivatives()
   {
      // Compute dkBasis
      int p, q, r, s, t;
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {

               // Loop over free indices s, t
               for (s = 0; s < D; ++s) {
                  for (t = 0; t < D; ++t) {
                     dkBasis_[p](q,r) 
                       -= kBasis_[q][s]*drBasis_[p](t,s)*kBasis_[t][r];
                  }
               }
               dkBasis_[p](q,r) /= 2.0*Constants::Pi;

            }
         }
      }

      // Compute drrBasis and dkkBasis 
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {
               for (s = 0; s < D; ++s) {
                  drrBasis_[p](q,r) += rBasis_[q][s]*drBasis_[p](r,s);
                  drrBasis_[p](q,r) += rBasis_[r][s]*drBasis_[p](q,s);
                  dkkBasis_[p](q,r) += kBasis_[q][s]*dkBasis_[p](r,s);
                  dkkBasis_[p](q,r) += kBasis_[r][s]*dkBasis_[p](q,s);
 
               }
            }
         }
      }

   }
   
   /*
   * Set all lattice parameters.
   */
   template <int D>
   void UnitCellBase<D>::setLattice()
   {
      initializeToZero();
      setBasis();
      computeDerivatives();
   }

} 
#endif 
