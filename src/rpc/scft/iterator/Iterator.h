#ifndef RPC_ITERATOR_H
#define RPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>  // base class template
#include <rpc/system/Types.h>           // base class template

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations in Rpc.
   *
   * \ingroup Rpc_Scft_Iterator_Module
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
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~Iterator() = default;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Iterator<1, Rpc::Types<1> >;
      extern template class Iterator<2, Rpc::Types<2> >;
      extern template class Iterator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Iterator<1>;
      extern template class Iterator<2>;
      extern template class Iterator<3>;
   }
} 
#endif
