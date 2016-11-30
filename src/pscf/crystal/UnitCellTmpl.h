#ifndef PSCF_UNIT_CELL_TMPL_H
#define PSCF_UNIT_CELL_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/FArray.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Template for a class representing a solvent species.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   class UnitCellTmpl 
   {
   public:

      /**
      * Constructor.
      */
      UnitCellTmpl()
      {}
   
      /**
      * Destructor.
      */
      ~UnitCellTmpl()
      {}
   
      /**
      * Get a Bravais basis vector i.
      */
      const RealVec<D>& rBasisVector(int i) const;
   
      /**
      * Get reciprocal basis vector i.
      */
      const RealVec<D>& kBasisVector(int i) const;
   
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
   const RealVec<D>& UnitCellTmpl<D>::rBasisVector(int i) const
   {  return rBasis_[i];  }

   /*
   * Get reciprocal basis vector i.
   */
   template <int D>
   inline
   const RealVec<D>& UnitCellTmpl<D>::kBasisVector(int i) const
   {  return kBasis_[i];  }

   /*
   * Get square magnitude of reciprocal basis vector.
   */
   template <int D>
   double UnitCellTmpl<D>::ksq(IntVec<D> const & k) const
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
