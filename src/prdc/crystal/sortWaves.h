#ifndef PRDC_SORT_WAVES_H
#define PRDC_SORT_WAVES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/Sort.h>
#include <vector>

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Prdc;

   /**
   * Generate sorted waves and identify bunches of equal magnitude.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void sortWaves(UnitCell<D> const & cell,
                  Mesh<D> const & dftMesh,
                  std::vector< Sort::Item<double> >& items,
                  std::vector< Sort::Slice >& bunches,
                  double epsilon);


   // Explicit instantiation declarations
   extern template 
   void sortWaves<1>(UnitCell<1> const &, Mesh<1> const &, 
                     std::vector< Sort::Item<double> >& , 
                     std::vector< Sort::Slice >& , double);

   extern template 
   void sortWaves<2>(UnitCell<2> const &, Mesh<2> const &, 
                     std::vector< Sort::Item<double> >& , 
                     std::vector< Sort::Slice >& , double);

   extern template 
   void sortWaves<3>(UnitCell<3> const &, Mesh<3> const &, 
                     std::vector< Sort::Item<double> >& , 
                     std::vector< Sort::Slice >& , double);

}
}
#endif
