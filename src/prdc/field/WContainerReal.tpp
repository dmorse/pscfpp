#ifndef PRDC_W_CONTAINER_REAL_TPP
#define PRDC_W_CONTAINER_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WContainerReal.h"
#include <prdc/field/fieldIoUtil.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <util/signal/Signal.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D, class RField, class FieldIo>
   WContainerReal<D,RField,FieldIo>::WContainerReal()
    : basis_(),
      rgrid_(),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      nMonomer_(0),
      signalPtr_(nullptr),
      fieldIoPtr_(nullptr),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false),
      hasData_(false),
      isSymmetric_(false)
   {
      signalPtr_ = new Signal<void>();
   }

   /*
   * Destructor.
   */
   template <int D, class RField, class FieldIo>
   WContainerReal<D,RField,FieldIo>::~WContainerReal()
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::setFieldIo(FieldIo const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D, class RField, class FieldIo>
   void WContainerReal<D,RField,FieldIo>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::allocateRGrid(IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);

      // If already allocated, deallocate.
      if (isAllocatedRGrid_) {
         deallocateRGrid();
      }

      // Store mesh dimensions
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }

      // Allocate arrays
      rgrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].allocate(meshDimensions);
      }
      isAllocatedRGrid_ = true;

   }

   /*
   * De-allocate memory for fields in r-grid format
   */
   template <int D, class RField, class FieldIo>
   void WContainerReal<D,RField,FieldIo>::deallocateRGrid()
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(nMonomer_ > 0);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].deallocate();
      }
      rgrid_.deallocate();
      meshDimensions_ = 0;
      meshSize_ = 0;
      isAllocatedRGrid_ = false;
   }

   /*
   * Allocate memory for fields in basis format.
   */
   template <int D, class RField, class FieldIo>
   void WContainerReal<D,RField,FieldIo>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(nBasis > 0);

      // If already allocated, deallocate.
      if (isAllocatedBasis_) {
         deallocateBasis();
      }

      // Allocate
      nBasis_ = nBasis;
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * De-allocate memory for fields in basis format.
   */
   template <int D, class RField, class FieldIo>
   void WContainerReal<D,RField,FieldIo>::deallocateBasis()
   {
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(nMonomer_ > 0);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].deallocate();
      }
      basis_.deallocate();
      nBasis_ = 0;
      isAllocatedBasis_ = false;
   }

   /*
   * Allocate memory for all fields.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::allocate(
                                        int nMonomer,
                                        int nBasis,
                                        IntVec<D> const & meshDimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(meshDimensions);
      allocateBasis(nBasis);
   }

   /*
   * Set new w-field values in symmetry-adapted basis format.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::setBasis(DArray< DArray<double> > const & fields)
   {
      UTIL_CHECK(fields.capacity() == nMonomer_);

      // Allocate fields if needed
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

      // Update system w fields (array basis_)
      for (int i = 0; i < nMonomer_; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> &  w = basis_[i];
         UTIL_CHECK(f.capacity() == nBasis_);
         UTIL_CHECK(w.capacity() == nBasis_);
         for (int j = 0; j < nBasis_; ++j) {
            w[j] = f[j];
         }
      }

      // Update system grid fields (array rgrid_)
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Set new field values, in r-grid format.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::setRGrid(DArray<RField> const & fields,
                                              bool isSymmetric)
   {
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         UTIL_CHECK(mesh.size() > 0);
         allocateRGrid(mesh.dimensions());
      }

      // Update rgrid_ fields
      UTIL_CHECK(fields.capacity() == nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         UTIL_CHECK(fields[i].capacity() == meshSize_);
         assignRField(rgrid_[i], fields[i]);
      }

      // If field isSymmetric, update basis fields
      if (isSymmetric) {
         if (!isAllocatedBasis_) {
            Basis<D> const & basis = fieldIo().basis();
            UTIL_CHECK(basis.isInitialized());
            allocateBasis(basis.nBasis());
         }
         fieldIo().convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ =  isSymmetric;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Read field component values from stream, in symmetrized basis format.
   *
   * This function also computes and stores the corresponding r-grid
   * representation. On return, hasData and isSymmetric are both true.
   */
   template <int D, class RField, class FieldIo>
   void 
   WContainerReal<D,RField,FieldIo>::readBasis(std::istream& in,
                                               UnitCell<D>& unitCell)
   {
      UTIL_CHECK(nMonomer_ > 0);

      // Read field file header
      int nMonomerIn;
      bool isSymmetricIn;
      fieldIo().readFieldHeader(in, nMonomerIn, unitCell, isSymmetricIn);
      // Note: FieldIo::readFieldHeader initializes basis if needed
      UTIL_CHECK(nMonomerIn == nMonomer_);
      UTIL_CHECK(isSymmetricIn);
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

      // Read data in basis form
      Prdc::readBasisData(in, basis_, unitCell, mesh, basis, nBasisIn);

      // Convert to r-grid form, to update rgrid_ array
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Read field component values from file, in symmetrized basis format.
   *
   * Calls readBasis(std::ifstream&, UnitCell<D>&) internally.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::readBasis(std::string filename,
                                               UnitCell<D>& unitCell)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readBasis(file, unitCell);
      file.close();
   }

   /*
   * Reads fields from an input stream in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the fields are known to be symmetric and so computes and stores 
   * the corresponding basis components. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   *
   * On return, hasData is true and the bool class member isSymmetric_ 
   * is set to the value of the isSymmetric function parameter.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::readRGrid(std::istream& in,
                                               UnitCell<D>& unitCell,
                                               bool isSymmetric)
   {
      // If necessary, allocate r-grid fields
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         allocateRGrid(mesh.dimensions());
      }

      // Read field file
      fieldIo().readFieldsRGrid(in, rgrid_, unitCell);

      // Optionally convert to symmetry-adapted basis form
      if (isSymmetric) {
         Basis<D> const & basis = fieldIo().basis();
         UTIL_CHECK(basis.isInitialized());
         if (!isAllocatedBasis_) {
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
   * Reads fields from a file in r-grid format, by filename.
   */
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::readRGrid(std::string filename,
                                               UnitCell<D>& unitCell,
                                               bool isSymmetric)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readRGrid(file, unitCell, isSymmetric);
      file.close();
   }

   /*
   * Symmetrize r-grid fields, convert to basis format.
   */
   template <int D, class RField, class FieldIo>
   void WContainerReal<D,RField,FieldIo>::symmetrize()
   {
      UTIL_CHECK(hasData_);
      fieldIo().convertRGridToBasis(rgrid_, basis_);
      fieldIo().convertBasisToRGrid(basis_, rgrid_);
      isSymmetric_ = true;

      // Notify observers of field modification
      signal().notify();
   }

   /*
   * Get a signal that is triggered by field modification.
   */
   template <int D, class RField, class FieldIo>
   Signal<void>& WContainerReal<D,RField,FieldIo>::signal()
   {
      UTIL_CHECK(signalPtr_);
      return *signalPtr_;
   }

   // Private virtual function

   /*
   * Assignment operation for r-grid fields (RField objects).
   *
   * Unimplemented virtual function - must be overridden by subclasses.
   */ 
   template <int D, class RField, class FieldIo>
   void
   WContainerReal<D,RField,FieldIo>::assignRField(RField & lhs,
                                                  RField const & rhs) const
   {  UTIL_THROW("Unimplemented function WContainerReal::assignRField");

}

} // namespace Prdc
} // namespace Pscf
#endif
