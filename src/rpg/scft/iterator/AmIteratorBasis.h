#ifndef RPG_AM_ITERATOR_BASIS_H
#define RPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.h>  // direct base class 
#include <rpg/system/Types.h>                  // direct base argument
#include <rpg/scft/iterator/Iterator.h>        // indirect base argument
#include <util/containers/DArray.h>            // indirect base argument

#if 0
namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;

   /**
   * Anderson mixing iterator with imposed space-group symmetry).
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AmIteratorBasis, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::AmIteratorBasis
   * \see \ref rp_AmIteratorBasis_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis : public Rp::AmIteratorBasis< D, Types<D> > 
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent system
      */
      AmIteratorBasis(System<D>& system);

   };

}
}
#endif

// Explicit instantiation declarations
namespace Pscf {
   extern template class AmIteratorTmpl<Rp::Iterator<1, Rpg::Types<1> >, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<2, Rpg::Types<2> >, DArray<double> >;
   extern template class AmIteratorTmpl<Rp::Iterator<3, Rpg::Types<3> >, DArray<double> >;
   namespace Rp {
      extern template class AmIteratorBasis<1, Rpg::Types<1> >;
      extern template class AmIteratorBasis<2, Rpg::Types<2> >;
      extern template class AmIteratorBasis<3, Rpg::Types<3> >;
   } 
}
#endif
