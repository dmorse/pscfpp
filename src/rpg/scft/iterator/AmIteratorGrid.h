#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorGrid.h>  // direct base class 
#include <rpg/system/Types.h>                 // direct base argument
#include <rpg/scft/iterator/Iterator.h>       // indirect base argument
#include <pscf/cuda/DeviceArray.h>            // indirect base argument
#include <pscf/cuda/cudaTypes.h>              // indirect base argument

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Anderson Mixing iterator on grid (no space-group symmetry).
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of base class template Rp::AmIteratorGrid, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::AmIteratorGrid
   * \see \ref rp_AmIteratorGrid_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page "AM Iteration Algorithm"
   * \ingroup Rpg_Scft_Iterator_Module
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

   extern template 
   class AmIteratorTmpl<Rp::Iterator<1, Rpg::Types<1> >, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl<Rp::Iterator<2, Rpg::Types<2> >, DeviceArray<cudaReal> >;
   extern template 
   class AmIteratorTmpl<Rp::Iterator<3, Rpg::Types<3> >, DeviceArray<cudaReal> >;

   namespace Rp {
      extern template class AmIteratorGrid<1, Rpg::Types<1> >;
      extern template class AmIteratorGrid<2, Rpg::Types<2> >;
      extern template class AmIteratorGrid<3, Rpg::Types<3> >;
   } 
   namespace Rpg {
      extern template class AmIteratorGrid<1>;
      extern template class AmIteratorGrid<2>;
      extern template class AmIteratorGrid<3>;
   } 
} 
#endif
