#ifndef PRDC_REPLICATE_UNIT_CELL_H
#define PRDC_REPLICATE_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // Template with a default parameter


namespace Pscf {
namespace Prdc {

   template <int D> class UnitCell;
   using namespace Util;
   using namespace Pscf;


   /**
   * Create a replicated UnitCell<D>.
   *
   * Element i of the IntVec<D> parameter named "replicas" contains the 
   * number of unit cell replicas along direction i. 
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of unit cell replicas in each direction
   * \param cellIn  original unit cell
   * \param cellOut  replicated unit cell
   */
   template <int D>
   void replicateUnitCell(IntVec<D> const & replicas,
                          UnitCell<D> const & cellIn,
                          UnitCell<D> & cellOut)
   {};

   /**
   * Create a replicated UnitCell<1>.
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of unit cell replicas in each direction
   * \param cellIn  original unit cell
   * \param cellOut  replicated unit cell
   */
   template<>
   void replicateUnitCell(IntVec<1> const & replicas,
                          UnitCell<1> const & cellIn,
                          UnitCell<1> & cellOut);

   /**
   * Create a replicated UnitCell<2>.
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of unit cell replicas in each direction
   * \param cellIn  original unit cell
   * \param cellOut  replicated unit cell
   */
   template<>
   void replicateUnitCell(IntVec<2> const & replicas,
                          UnitCell<2> const & cellIn,
                          UnitCell<2> & cellOut);

   /**
   * Create a replicated UnitCell<3>.
   *
   * Element i of the IntVec<D> parameter named "replicas" contains the 
   * number of unit cell replicas along direction i. 
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of unit cell replicas in each direction
   * \param cellIn  original unit cell 
   * \param cellOut  replicated unit cell (output)
   */
   template<>
   void replicateUnitCell(IntVec<3> const & replicas,
                          UnitCell<3> const & cellIn,
                          UnitCell<3> & cellOut);

} // namespace Prdc
} // namespace Pscf
#endif
