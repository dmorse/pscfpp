#ifndef PRDC_MASK_TMPL_TPP
#define PRDC_MASK_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskTmpl.h"
#include <prdc/field/fieldIoUtil.h> 
#include <prdc/crystal/Basis.h> 
#include <prdc/crystal/UnitCell.h> 
#include <pscf/mesh/Mesh.h> 
#include <util/signal/Signal.h>
#include <util/misc/FileMaster.h> 

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class FieldIo, class RField>
   MaskTmpl<D, FieldIo, RField>::MaskTmpl()
    : basis_(),
      rgrid_(),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      signalPtr_(nullptr),
      fieldIoPtr_(nullptr),
      isAllocatedBasis_(false),
      isAllocatedRGrid_(false),
      hasData_(false),
      isSymmetric_(false)
   {
      signalPtr_ = new Signal<void>();
   }

   /*
   * Destructor.
   */
   template <int D, class FieldIo, class RField>
   MaskTmpl<D,FieldIo,RField>::~MaskTmpl()
   {
      delete signalPtr_;
   }

   /*
   * Create an association with a FieldIo object.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::setFieldIo(FieldIo const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate memory for a field in basis format.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(!isAllocatedBasis_);

      // Set basis dimensions
      nBasis_ = nBasis;
  
      // Allocate field array, basis format 
      basis_.allocate(nBasis);
      isAllocatedBasis_ = true;
   }

   /*
   * Allocate memory for field in basis format.
   */
   template <int D, class FieldIo, class RField>
   void 
   MaskTmpl<D,FieldIo,RField>::allocateRGrid(IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(!isAllocatedRGrid_);

      // Set mesh dimensions
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate field array, rgrid format
      rgrid_.allocate(meshDimensions);
      isAllocatedRGrid_ = true;
   }

   /*
   * Set new field values, in basis form.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::setBasis(DArray<double> const & field)
   {
      // Allocate fields as needed
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         UTIL_CHECK(mesh.size() > 0);
         allocateRGrid(mesh.dimensions());
      }
      if (!isAllocatedBasis_) {
         Basis<D> const & basis = fieldIo().basis();
         UTIL_CHECK(basis.isInitialized());
         allocateBasis(basis.nBasis());
      }
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(isAllocatedBasis_);

      // Set components in basis form (array basis_)
      UTIL_CHECK(field.capacity() == nBasis_);
      for (int j = 0; j < nBasis_; ++j) {
         basis_[j] = field[j];
      }

      // Convert to r-grid form (update array rgrid_)
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Set new field values, in r-grid form.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::setRGrid(RField const & field,
                                             bool isSymmetric)
   {
      // Allocate rgrid_ field as needed
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         allocateRGrid(mesh.dimensions());
      }
      UTIL_CHECK(isAllocatedRGrid_);
   
      // Copy input field data to member variable rgrid_
      rgrid_ = field; // deep copy by assignment operator

      // Optionally convert to basis form
      if (isSymmetric) {
         if (!isAllocatedBasis_) {
            Basis<D> const & basis = fieldIo().basis();
            UTIL_CHECK(basis.isInitialized());
            allocateBasis(basis.nBasis());
         }
         UTIL_CHECK(isAllocatedBasis_);
         fieldIo().convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Read field from an input stream in basis format.
   *
   * This function also computes and stores the corresponding r-grid
   * representation. On return, hasData and isSymmetric are both true.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::readBasis(std::istream& in, 
                                              UnitCell<D>& unitCell)
   {
      // Read field file header
      int nMonomerIn;
      bool isSymmetricIn;
      fieldIo().readFieldHeader(in, nMonomerIn, unitCell, isSymmetricIn);
      UTIL_CHECK(1 == nMonomerIn);
      UTIL_CHECK(isSymmetricIn);
      UTIL_CHECK(fieldIo().basis().isInitialized());
      // Note: FieldIo::readFieldHeader will initialize basis if needed
      int nBasisIn = readNBasis(in);

      // Local references to mesh and basis
      Mesh<D> const & mesh = fieldIo().mesh();
      Basis<D> const & basis = fieldIo().basis();
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(basis.isInitialized());

      // If necessary, allocate fields
      if (!isAllocatedRGrid_) {
         allocateRGrid(mesh.dimensions());
      }
      if (!isAllocatedBasis_) {
         allocateBasis(basis.nBasis());
      }
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(isAllocatedBasis_);

      // Read field data in basis form (array basis_)
      Prdc::readBasisData(in, basis_, unitCell, mesh, basis, nBasisIn);

      // Convert r-grid form (array rgrid_)
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Read field components from a file basis format, by filename.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::readBasis(std::string filename, 
                                              UnitCell<D>& unitCell)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readBasis(file, unitCell);
      file.close();
   }

   /*
   * Reads field from an input stream in r-grid format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the field is known to be symmetric and so computes and stores
   * the corresponding basis format. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::readRGrid(std::istream& in, 
                                              UnitCell<D>& unitCell, 
                                              bool isSymmetric)
   {
      // If necessary, allocate rgrid_ fields
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         UTIL_CHECK(mesh.size() > 0);
         allocateRGrid(mesh.dimensions());
      }
      UTIL_CHECK(isAllocatedRGrid_);
  
      // Read field file in r-grid format (array rgrid_)
      fieldIo().readFieldRGrid(in, rgrid_, unitCell);

      // Optionally convert to basis form
      if (isSymmetric) {
         if (!isAllocatedBasis_) {
            Basis<D> const & basis = fieldIo().basis();
            UTIL_CHECK(basis.isInitialized());
            allocateBasis(basis.nBasis());
         }
         fieldIo().convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Read field from a file in r-grid format, by filename.
   */
   template <int D, class FieldIo, class RField>
   void MaskTmpl<D,FieldIo,RField>::readRGrid(std::string filename, 
                                              UnitCell<D>& unitCell, 
                                              bool isSymmetric)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readRGrid(file, unitCell, isSymmetric);
      file.close();
   }

   /*
   * Return volume fraction of the unit cell occupied by the 
   * polymers/solvents.
   */
   template <int D, class FieldIo, class RField>
   double MaskTmpl<D,FieldIo,RField>::phiTot() const
   {
      if (isSymmetric() && hasData()) {
         // Data in basis format is available
         return basis()[0];
      } else if (!hasData()) {
         // system does not have a mask
         return 1.0;
      } else { // Data is only available in r-grid format
         return rGridAverage();
      }
   }

   /*
   * Get a signal that is triggered by field modification.
   */
   template <int D, class FieldIo, class RField>
   Signal<void>& MaskTmpl<D,FieldIo,RField>::signal()
   {
      UTIL_CHECK(signalPtr_);
      return *signalPtr_;
   }

} // namespace Prdc
} // namespace Pscf
#endif
