#ifndef RPG_ITERATOR_H
#define RPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>  // base class template

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
   class Iterator : public Rp::Iterator<D, System<D> >
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
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~Iterator() = default;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Iterator<1, Rpg::System<1> >;
      extern template class Iterator<2, Rpg::System<2> >;
      extern template class Iterator<3, Rpg::System<3> >;
   }
   namespace Rpg {
      extern template class Iterator<1>;
      extern template class Iterator<2>;
      extern template class Iterator<3>;
   }
} 
#endif
