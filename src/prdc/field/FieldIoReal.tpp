#ifndef RPC_FIELD_IO_REAL_TPP
#define RPC_FIELD_IO_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIoReal.h"

#include <prdc/field/fieldIoUtil.h>
#include <prdc/crystal/fieldHeader.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/BFieldComparison.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class RFT, class KFT, class FFT>
   FieldIoReal<D,RFT,KFT,FFT>::FieldIoReal()
    : meshPtr_(nullptr),
      fftPtr_(nullptr),
      hasGroupPtr_(nullptr),
      groupNamePtr_(nullptr),
      groupPtr_(nullptr),
      basisPtr_(nullptr),
      fileMasterPtr_(nullptr),
      nMonomer_(0),
      isAllocatedBasis_(false),
      isAllocatedRGrid_(false),
      isAllocatedKGrid_(false)
   {}

   /*
   * Destructor.
   */
   template <int D, class RFT, class KFT, class FFT>
   FieldIoReal<D,RFT,KFT,FFT>::~FieldIoReal()
   {}

   // Initialization functions

   /*
   * Create associations with other members of a Domain<D> object.
   */
   template <int D, class RFT, class KFT, class FFT>
   void
   FieldIoReal<D,RFT,KFT,FFT>::associate(
                    Mesh<D> const & mesh,
                    FFT const & fft,
                    typename UnitCell<D>::LatticeSystem const & lattice,
                    bool const & hasGroup,
                    std::string const & groupName,
                    SpaceGroup<D> const & group,
                    Basis<D> & basis)
   {
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      latticePtr_ = &lattice;
      hasGroupPtr_ = &hasGroup;
      groupNamePtr_ = &groupName;
      groupPtr_ = &group;
      basisPtr_ = &basis;
   }

   /*
   * Create an association with a FileMaster.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::setFileMaster(
                              FileMaster const & fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Set nMonomer, the number of monomer types.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::setNMonomer(int nMonomer)
   {  nMonomer_ = nMonomer; }

   // Field File IO - Symmetry-Adapted Basis Format

   /*
   * Read an array of fields in basis format from an input stream.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldsBasis(
                              std::istream& in,
                              DArray< DArray<double> >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Precondition
      UTIL_CHECK(hasGroup());

      // Read header (checks compatibility with space group)
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(isSymmetric);
      UTIL_CHECK(basis().isInitialized());
      // Note: readFieldHeader can initialize basis if not done previously
      int nBasisIn = readNBasis(in);

      // Check allocation of fields container
      if (fields.isAllocated()) {
         int nMonomerFields, fieldCapacity;
         inspectArrays(fields, nMonomerFields, fieldCapacity);
         UTIL_CHECK(nMonomerFields == nMonomer);
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(basis().nBasis());
         }
      }

      // Read field data (components in basis)
      Prdc::readBasisData(in, fields, unitCell, mesh(), basis(), nBasisIn);
   }

   /*
   * Read a single field in basis format from an input stream
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldBasis(
                              std::istream& in,
                              DArray<double>& field,
                              UnitCell<D>& unitCell) const
   {
      // Local array container, of type required by readFieldsBasis
      DArray< DArray<double> > fields;

      // If single field is allocated, allocate local array fields
      if (field.isAllocated()) {
         fields.allocate(1);
         fields[0].allocate(field.capacity());
      }
      // Otherwise, pass unallocated fields array to readFieldsBasis

      // Read file containing a single field, allocate fields if needed.
      readFieldsBasis(in, fields, unitCell);

      // Check that only one field was read from file
      UTIL_CHECK(fields.capacity() == 1);

      // Copy data from local array fields to function parameter field
      field = fields[0];
   }

   /*
   * Write an array of fields in basis format to an output stream.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldsBasis(
                              std::ostream &out,
                              DArray< DArray<double> > const & fields,
                              UnitCell<D> const & unitCell) const
   {
      // Inspect fields to obtain nMonomer and fieldCapacity
      int nMonomer;
      int fieldCapacity;
      inspectArrays(fields, nMonomer, fieldCapacity);

      // Preconditions
      UTIL_CHECK(basis().isInitialized());
      UTIL_CHECK(fieldCapacity <= basis().nBasis());
      int nBasis = fieldCapacity;

      // Write header
      bool isSymmetric = true;
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeNBasis(out, nBasis);

      // Write data (field components)
      Prdc::writeBasisData(out, fields, basis());
   }

   /*
   * Write a single field in basis format to an output stream.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldBasis(
                              std::ostream& out,
                              DArray<double> const & field,
                              UnitCell<D> const & unitCell) const
   {
      // Create local array of type required by writeFieldsBasis
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());

      // Copy data from input parameter to local array
      fields[0] = field;

      writeFieldsBasis(out, fields, unitCell);
   }

   /*
   * File IO wrapper functions:
   *
   * These functions take a file name as an argument, and simply wrap
   * file open and close operations around a function of the same name
   * that takes an io stream argument. These functions can use the same
   * implementation in Rpc and Rpg.
   */

   /*
   * Open-close a file and read a set of fields in basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldsBasis(
                              std::string filename,
                              DArray<DArray<double> >& fields,
                              UnitCell<D>& unitCell) const
   {

      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsBasis(file, fields, unitCell);
      file.close();
   }

   /*
   * Open-close a file and read a single field in basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldBasis(
                              std::string filename,
                              DArray<double>& field,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldBasis(file, field, unitCell);
      file.close();
   }

   /*
   * Open-close a file, and write an array of fields in basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldsBasis(
                              std::string filename,
                              DArray<DArray<double> > const & fields,
                              UnitCell<D> const & unitCell) const
   {
       std::ofstream file;
       fileMaster().openOutputFile(filename, file);
       writeFieldsBasis(file, fields, unitCell);
       file.close();
   }

   /*
   * Open-close a file, and write a single field in basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldBasis(
                              std::string filename,
                              DArray<double> const & field,
                              UnitCell<D> const & unitCell) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldBasis(file, field, unitCell);
      file.close();
   }

   /*
   * Open-close a file, and write an array of fields in r-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   bool FieldIoReal<D,RFT,KFT,FFT>::readFieldsRGrid(
                              std::string filename,
                              DArray< RFT >& fields,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      bool isSymmetric;
      isSymmetric = readFieldsRGrid(file, fields, unitCell);
      file.close();
      return isSymmetric;
   }

   /*
   * Open-close a file, and read a single field in r-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   bool FieldIoReal<D,RFT,KFT,FFT>::readFieldRGrid(
                              std::string filename,
                              RFT & field,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      bool isSymmetric;
      isSymmetric = readFieldRGrid(file, field, unitCell);
      file.close();
      return isSymmetric;
   }

   /*
   * Open-close a file, and write an array of fields in r-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldsRGrid(
                              std::string filename,
                              DArray< RFT > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool writeHeader = true;
      bool writeMeshSize = true;
      writeFieldsRGrid(file, fields, unitCell,
                       writeHeader, isSymmetric, writeMeshSize);
      file.close();
   }

   /*
   * Open-close a file, and write a single field in r-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldRGrid(
                              std::string filename,
                              RFT const & field,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field, unitCell, isSymmetric);
      file.close();
   }

   /*
   * Open-close a file, and read an array of fields in k-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldsKGrid(
                              std::string filename,
                              DArray< KFT >& fields,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   /*
   * Open-close a file, and read an array of fields in k-grid format
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldsKGrid(
                              std::string filename,
                              DArray< KFT > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell, isSymmetric);
      file.close();
   }

   // Field Format Conversion Functions - Basis <-> KGrid

   /*
   * Convert array of fields from basis to k-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertBasisToKGrid(
                              DArray< DArray <double> > const & in,
                              DArray< KFT >& out) const
   {
      // Inspect input and output field containers
      int nMonomer, nMonomerOut, capacity;
      inspectArrays(in, nMonomer, capacity);
      IntVec<D> dimensions;
      inspectFields(out, nMonomerOut, dimensions);
      UTIL_CHECK(nMonomer == nMonomerOut);

      // Convert fields for all monomer types
      for (int i = 0; i < nMonomer; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   /*
   * Convert array of fields from k-grid format to basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertKGridToBasis(
                              DArray< KFT > const & in,
                              DArray< DArray <double> > & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      // Inspect input and output containers
      int nMonomer, nMonomerOut, capacity;
      IntVec<D> dimensions;
      inspectFields(in, nMonomer, dimensions);
      inspectArrays(out, nMonomerOut, capacity);
      UTIL_CHECK(nMonomer == nMonomerOut);

      // Convert fields for all monomer types
      for (int i = 0; i < nMonomer; ++i) {
         convertKGridToBasis(in[i], out[i], checkSymmetry, epsilon);
      }
   }

   /*
   * Convert a field file from k-grid to basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertKGridToBasis(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateKGrid();
      checkAllocateBasis(inFileName);
      UnitCell<D> tmpUnitCell;
      readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      convertKGridToBasis(tmpFieldsKGrid_, tmpFieldsBasis_);
      writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert a field file from basis to k-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertBasisToKGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateKGrid();
      checkAllocateBasis(inFileName);
      UnitCell<D> tmpUnitCell;
      readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      convertBasisToKGrid(tmpFieldsBasis_, tmpFieldsKGrid_);
      writeFieldsKGrid(outFileName, tmpFieldsKGrid_, tmpUnitCell);
   }

   // Field Format Conversion Functions - Basis <-> RGrid

   /*
   * Convert a single field from basis to r-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertBasisToRGrid(
                              DArray<double> const & in,
                              RFT& out) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      checkAllocateField(workDft_, mesh().dimensions());
      convertBasisToKGrid(in, workDft_);
      fft().inverseTransformUnsafe(workDft_, out);
   }

   /*
   * Convert an array of fields from basis to r-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertBasisToRGrid(
                              DArray< DArray <double> > const & in,
                              DArray< RFT >& out) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() == out.capacity());
      checkAllocateField(workDft_, mesh().dimensions());

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransformUnsafe(workDft_, out[i]);
      }
   }

   /*
   * Convert a single field from r-grid to basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToBasis(
                              RFT const & in,
                              DArray<double> & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      convertKGridToBasis(workDft_, out, checkSymmetry, epsilon);
   }

   /*
   * Convert an array of fields from r-grid to basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToBasis(
                              DArray< RFT > const & in,
                              DArray< DArray <double> > & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      // Inspect input and output containers
      int nMonomer, nMonomerOut, capacity;
      IntVec<D> dimensions;
      inspectFields(in, nMonomer, dimensions);
      UTIL_CHECK(mesh().dimensions() == dimensions);
      inspectArrays(out, nMonomerOut, capacity);
      UTIL_CHECK(nMonomer == nMonomerOut);
      checkAllocateField(workDft_, dimensions);

      // Convert RGrid -> KGrid -> Basis for each field
      for (int i = 0; i < nMonomer; ++i) {
         fft().forwardTransform(in[i], workDft_);
         convertKGridToBasis(workDft_, out[i], checkSymmetry, epsilon);
      }
   }

   /*
   * Convert a field file from basis to r-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertBasisToRGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateRGrid();
      checkAllocateBasis(inFileName);
      UnitCell<D> tmpUnitCell;

      readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      convertBasisToRGrid(tmpFieldsBasis_, tmpFieldsRGrid_);
      writeFieldsRGrid(outFileName, tmpFieldsRGrid_, tmpUnitCell);
   }

   /*
   * Convert a field file from r-grid to basis format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToBasis(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateRGrid();
      checkAllocateBasis(inFileName);
      UnitCell<D> tmpUnitCell;
      readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      convertRGridToBasis(tmpFieldsRGrid_, tmpFieldsBasis_);
      writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   // FFT conversion functions (KGrid <-> RGrid) [Same in Rpc and Rpg]

   /*
   * Apply inverse FFT to an array of k-grid fields, converting to r-grid.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertKGridToRGrid(
                              DArray< KFT > const & in,
                              DArray< RFT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   /*
   * Apply inverse FFT to a single k-grid field, converting to r-grid.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertKGridToRGrid(
                              KFT const & in, RFT& out) const
   {
      fft().inverseTransformSafe(in, out);
   }

   /*
   * Apply forward FFT to an array of r-grid fields, converting to k-grid.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToKGrid(
                              DArray< RFT > const & in,
                              DArray< KFT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Apply forward FFT to a single r-grid field, converting to k-grid.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToKGrid(
                              RFT const & in,
                              KFT& out) const
   {  fft().forwardTransform(in, out); }

   /*
   * Convert a field file from k-grid to r-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertKGridToRGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateKGrid();
      checkAllocateRGrid();
      UnitCell<D> tmpUnitCell;
      readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      for (int i = 0; i < nMonomer_; ++i) {
         fft().inverseTransformUnsafe(tmpFieldsKGrid_[i],
                                      tmpFieldsRGrid_[i]);
      }
      writeFieldsRGrid(outFileName, tmpFieldsRGrid_, tmpUnitCell);
   }

   /*
   * Convert a field file from r-grid to k-grid format.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::convertRGridToKGrid(
                                std::string const & inFileName,
                                std::string const & outFileName) const
   {
      checkAllocateRGrid();
      checkAllocateKGrid();
      UnitCell<D> tmpUnitCell;
      readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      convertRGridToKGrid(tmpFieldsRGrid_, tmpFieldsKGrid_);
      writeFieldsKGrid(outFileName, tmpFieldsKGrid_, tmpUnitCell);
   }

   // Field Inspection

   /*
   * Test if a single r-grid field has declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D, class RFT, class KFT, class FFT>
   bool FieldIoReal<D,RFT,KFT,FFT>::hasSymmetry(
                              RFT const & in,
                              double epsilon,
                              bool verbose) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_, epsilon, verbose);
   }

   /*
   * Check if r-grid fields have declared space group symmetry.
   */
   template <int D, class RFT, class KFT, class FFT>
   bool FieldIoReal<D,RFT,KFT,FFT>::hasSymmetry(
                                  std::string const & inFileName,
                                  double epsilon) const
   {
      checkAllocateRGrid();

      // Read fields
      UnitCell<D> tmpUnitCell;
      readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);

      // Check symmetry for all fields
      for (int i = 0; i < nMonomer_; ++i) {
         bool symmetric;
         symmetric = hasSymmetry(tmpFieldsRGrid_[i], epsilon);
         if (!symmetric) {
            return false;
         }
      }
      return true;
   }

   /*
   * Compare two fields in basis format, write report to Log file.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::compareFieldsBasis(
                        DArray< DArray<double> > const & field1,
                        DArray< DArray<double> > const & field2) const
   {
      BFieldComparison comparison(1);
      comparison.compare(field1, field2);

      Log::file() << "\n Basis expansion field comparison results"
                  << std::endl;
      Log::file() << "     Maximum Absolute Difference:   "
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   "
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::compareFieldsBasis(
                                std::string const & filename1,
                                std::string const & filename2) const
   {
      DArray< DArray<double> > fields1, fields2;
      UnitCell<D> tmpUnitCell;
      // Unallocated arrays will be allocated in readFieldsBasis
      readFieldsBasis(filename1, fields1, tmpUnitCell);
      readFieldsBasis(filename2, fields2, tmpUnitCell);
      compareFieldsBasis(fields1, fields2);
   }

   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::compareFieldsRGrid(
                                std::string const & filename1,
                                std::string const & filename2) const
   {
      DArray< RFT > fields1, fields2;
      UnitCell<D> tmpUnitCell;
      // Unallocated arrays will be allocated in readFieldsRGrid
      readFieldsRGrid(filename1, fields1, tmpUnitCell);
      readFieldsRGrid(filename2, fields2, tmpUnitCell);
      compareFieldsRGrid(fields1, fields2);
   }

   // Field Scaling

   /*
   * Multiply a single field in basis format by a constant factor.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::scaleFieldBasis(
                              DArray<double> & field,
                              double factor) const
   {
      UTIL_CHECK(field.isAllocated());
      int n = field.capacity();
      for (int i = 0; i < n; ++i) {
         field[i] *= factor;
      }
   }

   /*
   * Rescale an array of fields in basis format by a constant factor.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::scaleFieldsBasis(
                              DArray< DArray<double> >& fields,
                              double factor) const
   {
      int n = fields.capacity();
      UTIL_CHECK(n > 0);
      int m = fields[0].capacity();
      for (int i = 0; i < n; ++i) {
         UTIL_CHECK(fields[i].capacity() == m);
         scaleFieldBasis(fields[i], factor);
      }
   }

   /*
   * Rescale fields in files by a constant factor.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::scaleFieldsBasis(
                                std::string const & inFileName,
                                std::string const & outFileName,
                                double factor) const
   {
      checkAllocateBasis(inFileName);
      UnitCell<D> tmpUnitCell;
      readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      scaleFieldsBasis(tmpFieldsBasis_, factor);
      writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Rescale fields in r-grid format by a constant factor.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::scaleFieldsRGrid(
                              DArray< RFT > & fields,
                              double factor) const
   {
      int n = fields.capacity();
      for (int i = 0; i < n; ++i) {
         scaleFieldRGrid(fields[i], factor);
      }
   }

   /*
   * Rescale fields by a constant factor, read and write to file.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::scaleFieldsRGrid(
                                std::string const & inFileName,
                                std::string const & outFileName,
                                double factor) const
   {
      checkAllocateRGrid();
      UnitCell<D> tmpUnitCell;
      bool isSymmetric;
      isSymmetric = readFieldsRGrid(inFileName, tmpFieldsRGrid_,
                                    tmpUnitCell);
      scaleFieldsRGrid(tmpFieldsRGrid_, factor);
      writeFieldsRGrid(outFileName, tmpFieldsRGrid_, tmpUnitCell,
                       isSymmetric);
   }

   // Grid manipulation utilities

   /*
   * Replicate unit cell a specified number of times in each direction.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::replicateUnitCell(
                              std::string filename,
                              DArray<RFT> const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      replicateUnitCell(file, fields, unitCell, replicas);
      file.close();
   }

   /*
   * Replicate unit cell a specified number of times in each direction.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::replicateUnitCell(
                                std::string const & inFileName,
                                std::string const & outFileName,
                                IntVec<D> const & replicas) const
   {
      checkAllocateRGrid();
      UnitCell<D> tmpUnitCell;
      readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      replicateUnitCell(outFileName, tmpFieldsRGrid_, tmpUnitCell,
                        replicas);
   }

   /*
   * Expand the number of spatial dimensions of an RField.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::expandRGridDimension(
                              std::string filename,
                              DArray<RFT> const & fields,
                              UnitCell<D> const & unitCell, int d,
                              DArray<int> newGridDimensions) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      expandRGridDimension(file, fields, unitCell, d, newGridDimensions);
      file.close();
   }

   /*
   * Expand the number of spatial dimensions of an RField.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::expandRGridDimension(
                                std::string const & inFileName,
                                std::string const & outFileName,
                                int d,
                                DArray<int> newGridDimensions) const
   {
      checkAllocateRGrid();
      UnitCell<D> tmpUnitCell;
      readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      expandRGridDimension(outFileName, tmpFieldsRGrid_, tmpUnitCell,
                           d, newGridDimensions);
   }

   // File Header IO Utilities

   /*
   * Read common part of field header.
   *
   * Extracts number of monomers (i.e., number of fields) and unitCell
   * data (lattice type and parameters) from the file, and returns these
   * via non-const reference parameters nMonomer and unitCell.
   *
   * Also validates lattice type and groupName (if any), and constructs
   * a symmetry-adapted basis if there is a group but no basis.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::readFieldHeader(
                              std::istream& in,
                              int& nMonomer,
                              UnitCell<D>& unitCell,
                              bool & isSymmetric) const
   {
      // Preconditions
      UTIL_CHECK(latticePtr_);
      if (unitCell.lattice() == UnitCell<D>::Null) {
         UTIL_CHECK(unitCell.nParameter() == 0);
      } else {
         UTIL_CHECK(unitCell.nParameter() > 0);
         UTIL_CHECK(unitCell.lattice() == lattice());
      }

      // Read field header to set unitCell, groupNameIn, nMonomer
      int ver1, ver2;
      std::string groupNameIn;

      Pscf::Prdc::readFieldHeader(in, ver1, ver2, unitCell,
                                  groupNameIn, nMonomer);
      // Note: Function definition in prdc/crystal/fieldHeader.tpp

      // Checks of data from header
      UTIL_CHECK(ver1 == 1);
      //UTIL_CHECK(ver2 == 0);
      UTIL_CHECK(unitCell.isInitialized());
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.nParameter() > 0);

      // Validate or initialize lattice type
      if (lattice() != unitCell.lattice()) {
         Log::file() << std::endl
               << "Error - Mismatched lattice types "
               << "in FieldIo function readFieldHeader:\n"
               << "  FieldIo::lattice  :" << lattice() << "\n"
               << "  Unit cell lattice :" << unitCell.lattice()
               << "\n";
         UTIL_THROW("Mismatched lattice types");
      }

      // Check for presence of group name in the field file header
      isSymmetric = false;
      if (groupNameIn != "") {
         isSymmetric = true;
      }

      // Process group and basis (if any)
      if (hasGroup()) {

         // Check consistency of groupName values
         if (isSymmetric) {
            UTIL_CHECK(groupNamePtr_);
            if (groupNameIn != groupName()) {
               Log::file() << std::endl
                  << "Error - Mismatched group names in "
                  << "FieldIo member function readFieldHeader:\n"
                  << "  FieldIo::groupName :" << groupName() << "\n"
                  << "  Field file header  :" << groupNameIn << "\n";
               UTIL_THROW("Mismatched group names");
            }
         }

         // If there is a group but no basis, construct a basis
         UTIL_CHECK(basisPtr_);
         if (!basis().isInitialized()) {
            basisPtr_->makeBasis(mesh(), unitCell, group());
         }
         UTIL_CHECK(basis().isInitialized());

      } else {

         if (isSymmetric) {
            Log::file() << std::endl
               << "Warning: Group name found in a field file header"
               << "but no group declared in the parameter file.\n";
         }

      }

   }

   /*
   * Write a field file header.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::writeFieldHeader(
                              std::ostream &out,
                              int nMonomer,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      if (isSymmetric) {
         UTIL_CHECK(hasGroup());
         gName = groupName();
      }
      Pscf::Prdc::writeFieldHeader(out, v1, v2, unitCell,
                                   gName, nMonomer);
      // Note: This function is defined in prdc/crystal/fieldHeader.tpp
   }

   // Protected functions to check and allocate private workspace arrays

   /*
   * If necessary, allocate r-grid workspace.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::checkAllocateRGrid() const
   {
      if (isAllocatedRGrid_) return;

      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(mesh().size() > 0);
      IntVec<D> const & meshDimensions = mesh().dimensions();
      tmpFieldsRGrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         tmpFieldsRGrid_[i].allocate(meshDimensions);
      }
      isAllocatedRGrid_ = true;
   }

   /*
   * If necessary, allocate k-grid field workspace.
   */
   template <int D, class RFT, class KFT, class FFT>
   void FieldIoReal<D,RFT,KFT,FFT>::checkAllocateKGrid() const
   {
      if (isAllocatedKGrid_) return;

      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(mesh().size() > 0);
      IntVec<D> const & meshDimensions = mesh().dimensions();
      tmpFieldsKGrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         tmpFieldsKGrid_[i].allocate(meshDimensions);
      }
      isAllocatedKGrid_ = true;
   }

   /*
   * If necessary, allocate basis field workspace.
   */
   template <int D, class RFT, class KFT, class FFT>
   void
   FieldIoReal<D,RFT,KFT,FFT>::checkAllocateBasis(
                                   std::string const & inFileName) const
   {
      if (isAllocatedBasis_) return;

      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(hasGroup());
      if (!basis().isInitialized()) {
         // Peek at field header to initialize basis
         std::ifstream file;
         fileMaster().openInputFile(inFileName, file);
         int nMonomerIn;
         bool isSymmetricIn;
         UnitCell<D> tmpUnitCell;
         readFieldHeader(file, nMonomerIn, tmpUnitCell, isSymmetricIn);
         // Note: readFieldHeader can initialize basis if needed
         file.close();
      }
      UTIL_CHECK(basis().isInitialized());
      int nBasis = basis().nBasis();
      UTIL_CHECK(nBasis > 0);
      tmpFieldsBasis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         DArray<double>& field = tmpFieldsBasis_[i];
         field.allocate(nBasis);
         for (int j = 0; j < nBasis; ++j) {
            field[j] = 0.0;
         }
      }
      isAllocatedBasis_ = true;
   }

} // namespace Prdc
} // namespace Pscf
#endif
