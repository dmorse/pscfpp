#ifndef RPC_AM_ITERATOR_GRID_H
#define RPC_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.h>  // direct base class 
#include <rpc/system/Types.h>                 // direct base argument
#include <rpc/scft/iterator/Iterator.h>       // indirect base argument
#include <util/containers/DRArray.h>          // indirect base argument

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Anderson Mixing iterator on grid (no space-group symmetry).
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AmIteratorGrid, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::AmIteratorGrid
   * \see \ref rp_AmIteratorGrid_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   * \ingroup Rpc_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorGrid : public Rp::AmIteratorGrid<D, Types<D> > 
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent system
      */
      AmIteratorGrid(System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   extern template class AmIteratorTmpl<Rpc::Iterator<1>, DRArray<double> >;
   extern template class AmIteratorTmpl<Rpc::Iterator<2>, DRArray<double> >;
   extern template class AmIteratorTmpl<Rpc::Iterator<3>, DRArray<double> >;
   namespace Rp {
      extern template class AmIteratorGrid<1, Rpc::Types<1> >;
      extern template class AmIteratorGrid<2, Rpc::Types<2> >;
      extern template class AmIteratorGrid<3, Rpc::Types<3> >;
   } 
   namespace Rpc {
      extern template class AmIteratorGrid<1>;
      extern template class AmIteratorGrid<2>;
      extern template class AmIteratorGrid<3>;
   } 
}
#endif
