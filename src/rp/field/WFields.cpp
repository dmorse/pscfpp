/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <prdc/field/cpp/RField.h>
#include <util/containers/DArray.h>

#include <rp/field/WFields.tpp> 
//#include <rp/field/WFields_c.h>     // class specialization

// Explicit instantiation definitions
namespace Pscf {
namespace Rp {

   #if 0
   /*
   * Set new w-field values, using unfolded array of r-grid fields.
   */
   template <int D>
   void WFields<D,CPT>::setRGrid(DeviceArray<cudaReal,CPT>& fields)
   {
      // Create DArray tmp with RField<D,CPT> elements
      DArray< RField<D,CPT> > tmp;
      const int nMonomer = Base::nMonomer();
      tmp.allocate(nMonomer);

      // Associate each RField<D,CPT> with a slice of the unfolded array
      IntVec<D> const & meshDimensions = Base::meshDimensions();
      const int meshSize = Base::meshSize();
      for (int i = 0; i < nMonomer; i++) {
         tmp[i].associate(fields, i*meshSize, meshDimensions);
      }

      // Use tmp array to set w-fields for all monomer types
      bool isSymmetric = false;
      Base::setRGrid(tmp, isSymmetric);
   }

   // Explicit instantiation definitions - base class
   template class WFieldsBase<1,CPT>;
   template class WFieldsBase<2,CPT>;
   template class WFieldsBase<3,CPT>;
   #endif

   // Explicit instantiation definitions - this subclass
   template class WFields<1,CPT>;
   template class WFields<2,CPT>;
   template class WFields<3,CPT>;

}
}
