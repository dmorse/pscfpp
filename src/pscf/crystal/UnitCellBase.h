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
      UnitCellBase()
      {}
   
      /**
      * Destructor.
      */
      ~UnitCellBase()
      {}
   
      /**
      * Get a Bravais basis vector i.
      */
      const RealVec<D>& rBasis(int i) const;
   
      /**
      * Get reciprocal basis vector i.
      */
      const RealVec<D>& kBasis(int i) const;

      /** 
      * Get the parameters of unit cell.
      */  
      FArray<double, 6> params();

      /** 
      * Set the parameters of unit cell.
      */  
      void SetParams(double val, int m);

      /**
      * Get the jth component of ith direction of derivative of rBasis basis vector with respect to k.
      */
      const double drBasis(int k, int i, int j) const;

      /**
      * Get the jth component of ith direction of derivative of kBasis basis vector with respect to k.
      */
      const double dkBasis(int k, int i, int j) const;
 
      /**
      * Get the derivative of ki and kj with respect to parameter k.
      */
      const double dkkBasis(int k, int i, int j) const;

      /** 
      * Get the derivative of ri and rj with respect to parameter k.
      */  
      const double drrBasis(int k, int i, int j) const;
  
      /**
      * Get the number of Parameters in the unit cell.
      */
      const int nParams() const;

      /**
      * Get square magnitude of reciprocal basis vector.
      */
      virtual double ksq(IntVec<D> const & k) const;
   
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
      * Array of derivatives of rBasis
      */ 
      FArray<FMatrix<double, D, D>, 6> drBasis_;
      
      /**
      * Array of derivatives of kBasis
      */
      FArray<FMatrix<double, D, D>, 6> dkBasis_;

      /**
      * Array of derivatives of a_i.a_j
      */
      FArray<FMatrix<double, D, D>, 6> drrBasis_;

      /**
      * Array of derivatives of b_i.b_j
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

   };

   /*
   * Get a Bravais basis vector i.
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
   * Get the jth component of ith direction of derivative of rBasis basis vector with respect to k.
   */
   template <int D>
   inline
   const double UnitCellBase<D>::drBasis(int k, int i, int j) const
   {  return drBasis_[k](i,j);  }

   /*
   * Get the jth component of ith direction of derivative of kBasis basis vector with respect to k.
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
   * Get the number of Parameters in the unit cell.
   */
   template <int D>
   inline
   const int UnitCellBase<D>::nParams() const
   {  return nParameter_;  }

  /*
   * Get the parameters in the unit cell.
   */
   template <int D>
   inline
   FArray<double, 6> UnitCellBase<D>::params()
   {  return parameters_;  }

  /*
   * Set the parameters in the unit cell.
   */
   template <int D>
   inline
   void UnitCellBase<D>::SetParams(double val, int m)
   {  parameters_ [m] = val;}

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
   
} 
#endif 
