#ifndef PRDC_REPLICATE_UNIT_CELL_H
#define PRDC_REPLICATE_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // Template with a default parameter

namespace Pscf {
namespace Prdc {

   template <int D> class UnitCell;
   using namespace Util;
   using namespace Pscf;


   /**
   * Create a replicated UnitCell<D> (base template).
   *
   * Element i of the IntVec<D> parameter "replicas" contains the number
   * of unit cell replicas along direction i, for i=0, ..., D-1.
   *
   * Explicit specializations of this function are used to implement the
   * REPLICATE_UNIT_CELL command of pscf_pc and pscf_pg, which replicates 
   * the unit cell and the fields defined within it a specified number of 
   * times in each of D directions. Explicit specializations of this 
   * function for D=1, 2, and 3 create a UnitCell<D> object appropriate 
   * for the replicated system. It is used within the replicateUnitCell 
   * function defined in file src/prdc/fields/fieldIoUtil.tpp, which 
   * writes replicated fields to an output file. 
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of replicas in each direction (in)
   * \param cellIn  original unit cell (in)
   * \param cellOut  replicated unit cell (out)
   */
   template <int D>
   void replicateUnitCell(IntVec<D> const & replicas,
                          UnitCell<D> const & cellIn,
                          UnitCell<D> & cellOut)
   {};

   // Explicit specializations of replicateUnitCell<D>

   /**
   * Create a replicated UnitCell<1>.
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of replicas in each direction (in)
   * \param cellIn  original unit cell (in)
   * \param cellOut  replicated unit cell (out)
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
   * \param replicas  number of replicas in each direction (in)
   * \param cellIn  original unit cell (in)
   * \param cellOut  replicated unit cell (out)
   */
   template<>
   void replicateUnitCell(IntVec<2> const & replicas,
                          UnitCell<2> const & cellIn,
                          UnitCell<2> & cellOut);

   /**
   * Create a replicated UnitCell<3>.
   *
   * \ingroup Prdc_Crystal_Module
   * 
   * \param replicas  number of replicas in each direction (in)
   * \param cellIn  original unit cell (in)
   * \param cellOut  replicated unit cell (out)
   */
   template<>
   void replicateUnitCell(IntVec<3> const & replicas,
                          UnitCell<3> const & cellIn,
                          UnitCell<3> & cellOut);

} // namespace Prdc
} // namespace Pscf
#endif
