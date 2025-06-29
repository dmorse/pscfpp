#ifndef PRDC_W_FIELDS_REAL_TPP
#define PRDC_W_FIELDS_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldsReal.h"
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
   template <int D, class RFT, class FIT>
   WFieldsReal<D,RFT,FIT>::WFieldsReal()
    : basis_(),
      rgrid_(),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      nMonomer_(0),
      readUnitCellPtr_(nullptr),
      writeUnitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      signalPtr_(nullptr),
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
   template <int D, class RFT, class FIT>
   WFieldsReal<D,RFT,FIT>::~WFieldsReal()
   {
      delete signalPtr_;
   }

   /*
   * Create an association with a FIT object.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::setFieldIo(FIT const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D, class RFT, class FIT>
   void WFieldsReal<D,RFT,FIT>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Set the unit cell that is modified by reading a field file.
   */
   template <int D, class RFT, class FIT>
   void WFieldsReal<D,RFT,FIT>::setReadUnitCell(UnitCell<D>& cell)
   {
      UTIL_CHECK(!readUnitCellPtr_);
      readUnitCellPtr_ = &cell;
   }

   /*
   * Set the unit cell that whose parameters are written to a field header.
   */
   template <int D, class RFT, class FIT>
   void WFieldsReal<D,RFT,FIT>::setWriteUnitCell(UnitCell<D> const & cell)
   {
      UTIL_CHECK(!writeUnitCellPtr_);
      writeUnitCellPtr_ = &cell;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::allocateRGrid(IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!hasData_);
      UTIL_CHECK(!isAllocatedRGrid_);

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
   * Allocate memory for fields in basis format.
   */
   template <int D, class RFT, class FIT>
   void WFieldsReal<D,RFT,FIT>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(nBasis > 0);
      UTIL_CHECK(!isAllocatedBasis_);
      UTIL_CHECK(!hasData_);

      nBasis_ = nBasis;
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * Allocate memory for all fields.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::allocate(int nMonomer,
                                       int nBasis,
                                       IntVec<D> const & meshDimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(meshDimensions);
      allocateBasis(nBasis);
   }

   // Field Modification Functions

   /*
   * Set new field values, in basis form.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::setBasis(DArray< DArray<double> > const & fields)
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

      // Set components in basis form (array basis_)
      for (int i = 0; i < nMonomer_; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> &  w = basis_[i];
         UTIL_CHECK(f.capacity() == nBasis_);
         UTIL_CHECK(w.capacity() == nBasis_);
         for (int j = 0; j < nBasis_; ++j) {
            w[j] = f[j];
         }
      }

      // Convert to r-grid form (update array rgrid_)
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify signal observers of field modification
      signal().notify();
   }

   /*
   * Set new field values, in r-grid form.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::setRGrid(DArray<RFT> const & fields,
                                       bool isSymmetric)
   {
      // Allocate r-grid fields as needed
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         UTIL_CHECK(mesh.size() > 0);
         allocateRGrid(mesh.dimensions());
      }
      UTIL_CHECK(isAllocatedRGrid_);

      // Update rgrid_ fields
      UTIL_CHECK(fields.capacity() == nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         UTIL_CHECK(fields[i].capacity() == meshSize_);
         assignRField(rgrid_[i], fields[i]);
      }

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
      isSymmetric_ =  isSymmetric;

      // Notify signal observers of field modification
      signal().notify();
   }

   /*
   * Read fields from an input stream in basis format.
   *
   * This function also computes and stores the corresponding r-grid
   * representation. On return, hasData and isSymmetric are both true.
   */
   template <int D, class RFT, class FIT>
   void 
   WFieldsReal<D,RFT,FIT>::readBasis(std::istream& in)
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(readUnitCellPtr_);

      // Read field file header
      int nMonomerIn;
      bool isSymmetricIn;
      fieldIo().readFieldHeader(in, nMonomerIn, *readUnitCellPtr_, 
                                isSymmetricIn);
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

      // Read data in basis form (array basis_)
      Prdc::readBasisData(in, basis_, 
                          *readUnitCellPtr_, mesh, basis, nBasisIn);

      // Convert to r-grid form (array rgrid_)
      fieldIo().convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;

      // Notify signal observers of field modification
      signal().notify();
   }

   /*
   * Read fields from a file in basis format, by filename.
   *
   * Calls readBasis(std::ifstream&) internally.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::readBasis(std::string filename)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readBasis(file);
      file.close();
   }

   /*
   * Read fields from an input stream in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the fields are known to be symmetric and so computes and stores 
   * the corresponding basis components. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   *
   * On return, hasData is true and the bool class member isSymmetric_ 
   * is set to the value of the isSymmetric function parameter.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::readRGrid(std::istream& in, 
                                        bool isSymmetric)
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(readUnitCellPtr_);
      
      // If necessary, allocate r-grid fields
      if (!isAllocatedRGrid_) {
         Mesh<D> const & mesh = fieldIo().mesh();
         allocateRGrid(mesh.dimensions());
      }

      // Read field file in r-grid format (array rgrid_)
      fieldIo().readFieldsRGrid(in, rgrid_, *readUnitCellPtr_);

      // Optionally convert to basis form
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

      // Notify signal observers of field modification
      signal().notify();
   }

   /*
   * Read fields from a file in r-grid format, by filename.
   */
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::readRGrid(std::string filename,
                                        bool isSymmetric)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readRGrid(file, isSymmetric);
      file.close();
   }

   /*
   * Symmetrize r-grid fields, convert to basis format.
   */
   template <int D, class RFT, class FIT>
   void WFieldsReal<D,RFT,FIT>::symmetrize()
   {
      UTIL_CHECK(hasData_);
      fieldIo().convertRGridToBasis(rgrid_, basis_);
      fieldIo().convertBasisToRGrid(basis_, rgrid_);
      isSymmetric_ = true;

      // Notify signal observers of field modification
      signal().notify();
   }

   /*
   * Get the Signal<void> that is triggered by field modification.
   */
   template <int D, class RFT, class FIT>
   Signal<void>& WFieldsReal<D,RFT,FIT>::signal()
   {
      UTIL_CHECK(signalPtr_);
      return *signalPtr_;
   }

   // Private virtual function

   /*
   * Assignment operation for r-grid fields (RFT objects).
   *
   * Unimplemented virtual function - must be overridden by subclasses.
   */ 
   template <int D, class RFT, class FIT>
   void
   WFieldsReal<D,RFT,FIT>::assignRField(RFT & lhs, RFT const & rhs) const
   {  UTIL_THROW("Unimplemented function WFieldsReal::assignRField");

}

} // namespace Prdc
} // namespace Pscf
#endif
