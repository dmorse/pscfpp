/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "sortWaves.h"
#include <prdc/cpu/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
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
                  IntVec<D> const & meshDimensions,
                  std::vector< Sort::Item<double> >& items,
                  std::vector< Sort::Bunch >& bunches,
                  double epsilon,
                  bool isRealField)
   {
      // Compute dimensions of Fourier space (DFT) k-grid
      IntVec<D> kMeshDimensions;
      int kMeshSize;
      if (isRealField) {
         Cpu::FFT<D>::computeKMesh(meshDimensions,
                                   kMeshDimensions, kMeshSize);
      } else {
         kMeshSize = 1;
         for (int i = 0; i < D; ++i) {
            kMeshDimensions[i] = meshDimensions[i];
            kMeshSize *= meshDimensions[i];
         }
      }

      // Generated unsorted array of items, one per wavevector
      items.clear();
      items.reserve(kMeshSize);
      Sort::Item<double> item;
      IntVec<D> v;
      MeshIterator<D> itr(kMeshDimensions);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         item.id = itr.rank();
         v = itr.position();
         v = shiftToMinimum(v, meshDimensions, cell);
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
   void sortWaves<1>(UnitCell<1> const &, IntVec<1> const &,
                     std::vector< Sort::Item<double> >&,
                     std::vector< Sort::Bunch >&, double, bool);
   template
   void sortWaves<2>(UnitCell<2> const &, IntVec<2> const &,
                     std::vector< Sort::Item<double> >&,
                     std::vector< Sort::Bunch >&, double, bool);
   template
   void sortWaves<3>(UnitCell<3> const &, IntVec<3> const &,
                     std::vector< Sort::Item<double> >&,
                     std::vector< Sort::Bunch >&, double, bool);

}
}
