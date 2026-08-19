#ifndef PRDC_SORT_WAVES_H
#define PRDC_SORT_WAVES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/Sort.h>
#include <pscf/math/IntVec.h>
#include <util/containers/GArray.h>
#include <vector>

// Forward declarations
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Prdc;

   /**
   * Generate sorted waves and identify "bunches" of equal magnitude.
   *
   * This function generates all wavevectors in the k-grid mesh associated
   * with an r-grid mesh of specified dimensions, sorts them in order of
   * increasing vector norm and identifies all "bunches" of wavevectors of 
   * equal norm. The function can generate wavevectors associated with the
   * k-grid mesh used for either complex or real-valued fields, depending
   * on the value of the bool parameter isRealField. 
   *
   * Upon return:
   *
   * Element i of array "items" is a Sort::Item<double> object for which 
   * items[i].value is the squared norm of an associated wavevector and
   * items[i].id is the rank of that wavevector within the k-grid mesh,
   * while these items are now sorted in order of increasing wavevector
   * norm.
   *
   * Element i of array "bunches" is a Sort::Bunch objects for which 
   * elements of the items array with indices between bunches[i][0]
   * and (bunches[i][1] - 1) have equal norms. Every bunch has one or
   * more associated wavevectors, and every wavevector belongs to a
   * unique bunch. 
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param cell  crystallographic unit cell (in)
   * \param meshDimensions  dimensions of real space mesh (in)
   * \param items  sorted array of items associated with wavevectors (out)
   * \param bunches  slices of items with equal squared wavenumbers (out)
   * \param epsilon  tolerance used to test for equality (in)
   * \param isRealField  Use k-grid for DFT of a real field? (in)
   */
   template <int D>
   void sortWaves(UnitCell<D> const & cell,
                  IntVec<D> const & meshDimensions,
                  std::vector< Sort::Item<double> >& items,
                  GArray< Sort::Bunch >& bunches,
                  double epsilon,
                  bool isRealField = true);

   // Explicit instantiation declarations
   extern template 
   void sortWaves<1>(UnitCell<1> const &, IntVec<1> const &, 
                     std::vector< Sort::Item<double> >& , 
                     GArray< Sort::Bunch >& , double, bool);

   extern template 
   void sortWaves<2>(UnitCell<2> const &, IntVec<2> const &, 
                     std::vector< Sort::Item<double> >& , 
                     GArray< Sort::Bunch >& , double, bool);

   extern template 
   void sortWaves<3>(UnitCell<3> const &, IntVec<3> const &, 
                     std::vector< Sort::Item<double> >& , 
                     GArray< Sort::Bunch >& , double, bool);

}
}
#endif
