#ifndef RPG_ITERATOR_H
#define RPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>  // base class template
#include <rpg/system/Types.h>           // base class template argument

#if 0
namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations in Rpg.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class Iterator : public Rp::Iterator<D, Types<D> >
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Iterator() = default;

   };

} // namespace Rpg
} // namespace Pscf
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Iterator<1, Rpg::Types<1> >;
      extern template class Iterator<2, Rpg::Types<2> >;
      extern template class Iterator<3, Rpg::Types<3> >;
   }
} 
#endif
