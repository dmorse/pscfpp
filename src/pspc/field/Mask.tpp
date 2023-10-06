#ifndef PSPC_MASK_TPP
#define PSPC_MASK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <pspc/field/FieldIo.h>
#include <prdc/cpu/RFieldDft.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Mask<D>::Mask()
    : basis_(),
      rgrid_(),
      fieldIoPtr_(0),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      isAllocated_(false),
      hasData_(false),
      isSymmetric_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   Mask<D>::~Mask()
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D>
   void Mask<D>::setFieldIo(FieldIo<D> const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate memory for field.
   */
   template <int D>
   void Mask<D>::allocate(int nBasis, IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(!isAllocated_);

      // Set mesh and basis dimensions
      nBasis_ = nBasis;
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate field arrays 
      basis_.allocate(nBasis);
      rgrid_.allocate(meshDimensions);
      isAllocated_ = true;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void Mask<D>::setBasis(DArray<double> const & field)
   {
      UTIL_CHECK(field.capacity() == nBasis_);
      for (int j = 0; j < nBasis_; ++j) {
         basis_[j] = field[j];
      }
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);
      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Set new field values, using r-grid field as inputs.
   */
   template <int D>
   void Mask<D>::setRGrid(RField<D> const & field,
                          bool isSymmetric)
   {
      for (int j = 0; j < meshSize_; ++j) {
         rgrid_[j] = field[j];
      }
      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }
      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

   /*
   * Read field from input stream, in symmetrized Fourier format.
   *
   * This function also computes and stores the corresponding
   * r-grid representation. On return, hasData and isSymmetric
   * are both true.
   */
   template <int D>
   void Mask<D>::readBasis(std::istream& in, UnitCell<D>& unitCell)
   {
      fieldIoPtr_->readFieldBasis(in, basis_, unitCell);

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Read field from file, in symmetrized Fourier format.
   *
   * This function also computes and stores the corresponding
   * r-grid representation. On return, hasData and isSymmetric
   * are both true.
   */
   template <int D>
   void Mask<D>::readBasis(std::string filename, UnitCell<D>& unitCell)
   {
      fieldIoPtr_->readFieldBasis(filename, basis_, unitCell);

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Reads field from an input stream in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the field is known to be symmetric and so computes and stores
   * the corresponding basis format. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   * 
   * On return, hasData is true and the persistent isSymmetric flag 
   * defined by the class is set to the value of the isSymmetric 
   * input parameter.
   */
   template <int D>
   void Mask<D>::readRGrid(std::istream& in, UnitCell<D>& unitCell, 
                           bool isSymmetric)
   {
      fieldIoPtr_->readFieldRGrid(in, rgrid_, unitCell);

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

   /*
   * Reads field from a file in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the field is known to be symmetric and so computes and stores
   * the corresponding basis format. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   * 
   * On return, hasData is true and the persistent isSymmetric flag 
   * defined by the class is set to the value of the isSymmetric 
   * input parameter.
   */
   template <int D>
   void Mask<D>::readRGrid(std::string filename, UnitCell<D>& unitCell, 
                           bool isSymmetric)
   {
      fieldIoPtr_->readFieldRGrid(filename, rgrid_, unitCell);

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

   /*
   * Return volume fraction of the unit cell occupied by the 
   * polymers/solvents.
   */
   template <int D>
   double Mask<D>::phiTot() const
   {
      if (isSymmetric() && hasData()) {
         // Data in basis format is available
         return basis()[0];
      } else if (!hasData()) {
         // system does not have a mask
         return 1.0;
      } else { // Data is only available in r-grid format
         RFieldDft<D> kgrid;
         kgrid.allocate(rgrid_.meshDimensions());
         fieldIoPtr_->convertRGridToKGrid(rgrid_, kgrid);
         UTIL_CHECK(kgrid[0][1] < 1e-8);
         return kgrid[0][0];
      }
   }

} // namespace Pspc
} // namespace Pscf
#endif
