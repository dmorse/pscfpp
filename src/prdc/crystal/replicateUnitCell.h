#ifndef PRDC_REPLICATE_UNIT_CELL_H
#define PRDC_REPLICATE_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * Upon return, parameter "cellOut" is a UnitCell<D> with a lattice
   * type and lattice parameters appropriate to represent a unit cell
   * created by replicating "cellIn" the specified number of times in 
   * each direction. This function attempts to identify symmetries that
   * are preserved by the choice of values in the "replicas" vector, and
   * will generally try to use the highest symmetry unit cell type that
   * is appropriate if more than one choice is possible.  Thus, if a
   * 3D cubic unit cell cellIn is replicated using a "replicas" vector 
   * with three equal elements, cellOut will be a cubic cell, whereas 
   * if a cubic cell replicated using a "replicas" vector in which all 
   * three elements are unequal, "cellOut" will be orthorhombic. 
   *
   * Some types of conversion of unit cells with high symmetry lattice
   * types that would result in lower-symmetry lattice type have been 
   * prohibited if the appropriate conversion of lattice parameters has
   * not yet been worked out. For example, for now, users may not
   * create a triclinic unit cell by replicating a rhombohedral lattice
   * using a replicas vector with unequal values for different elements.
   * An Exception is thrown and an appropriate error message is written 
   * to standard output if a user attempts any such prohibited type of
   * replication.
   *
   * This base template is declared but not defined. Explicit 
   * specializations are defined for D=1, 2, and 3.  
   * 
   * Explicit specializations of this function template are called 
   * within the replicateUnitCell function template defined in file 
   * src/prdc/fields/fieldIoUtil.tpp. This function takes an array of
   * fields as an input and writes a replicated version of those fields
   * to an r-grid field file. The header of the resulting file contains a 
   * description of the unit cell obtained here. That field replication 
   * function is used to implement the REPLICATE_UNIT_CELL command of
   * the pscf_pc and pscf_pg programs.
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
                          UnitCell<D> & cellOut);

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
