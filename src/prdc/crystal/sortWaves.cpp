/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "sortWaves.h"
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Sorted waves and identify bunches of equal magnitude.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void sortWaves(UnitCell<D> const & cell,
                  Mesh<D> const & dftMesh,
                  std::vector< Sort::Item<double> >& items,
                  std::vector< Sort::Slice >& bunches,
                  double epsilon)
   {
      IntVec<D> dftMeshDimensions = dftMesh.dimensions();
      int dftMeshSize = dftMesh.size();

      // Generated unsorted array of items, one per wavevector
      items.clear();
      items.reserve(dftMeshSize);
      Sort::Item<double> item;
      IntVec<D> v;
      MeshIterator<D> itr(dftMeshDimensions);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         item.id = itr.rank();
         v = itr.position();
         v = shiftToMinimum(v, dftMeshDimensions, cell);
         item.value = cell.ksq(v);
         items.push_back(item);
      }

      // Sort items
      Sort::sort(items);
   
      // Identify bunches (array slices) of waves of equal magnitude 
      bunches.clear(); 
      Sort::findBunches(items, bunches, epsilon); 
   }

   // Explicit instantiation definitions
   template 
   void sortWaves<1>(UnitCell<1> const &, Mesh<1> const &, 
                     std::vector< Sort::Item<double> >&, 
                     std::vector< Sort::Slice >&, double);
   template 
   void sortWaves<2>(UnitCell<2> const &, Mesh<2> const &, 
                     std::vector< Sort::Item<double> >&, 
                     std::vector< Sort::Slice >&, double);
   template 
   void sortWaves<3>(UnitCell<3> const &, Mesh<3> const &, 
                     std::vector< Sort::Item<double> >&, 
                     std::vector< Sort::Slice >&, double);

}
}
