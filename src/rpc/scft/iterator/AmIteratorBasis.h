#ifndef RPC_AM_ITERATOR_BASIS_H
#define RPC_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.h>  // direct base class 
#include <rpc/system/Types.h>                  // direct base argument
#include <rpc/scft/iterator/Iterator.h>        // indirect base argument
#include <util/containers/DArray.h>            // indirect base argument

namespace Pscf {
namespace Rpc {

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
   * \ingroup Rpc_Scft_Iterator_Module
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

// Explicit instantiation declarations
namespace Pscf {
   extern template class AmIteratorTmpl<Rpc::Iterator<1>, DArray<double> >;
   extern template class AmIteratorTmpl<Rpc::Iterator<2>, DArray<double> >;
   extern template class AmIteratorTmpl<Rpc::Iterator<3>, DArray<double> >;
   namespace Rp {
      extern template class AmIteratorBasis<1, Rpc::Types<1> >;
      extern template class AmIteratorBasis<2, Rpc::Types<2> >;
      extern template class AmIteratorBasis<3, Rpc::Types<3> >;
   } 
   namespace Rpc {
      extern template class AmIteratorBasis<1>;
      extern template class AmIteratorBasis<2>;
      extern template class AmIteratorBasis<3>;
   } 
}
#endif
