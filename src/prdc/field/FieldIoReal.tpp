#ifndef RPC_FIELD_IO_REAL_TPP
#define RPC_FIELD_IO_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   template <int D, class RFRT, class RFKT, class FFTT>
   FieldIoReal<D,RFRT,RFKT,FFTT>::FieldIoReal()
    : meshPtr_(0),
      fftPtr_(0),
      hasGroupPtr_(0),
      groupNamePtr_(0),
      groupPtr_(0),
      basisPtr_(0),
      fileMasterPtr_()
   {}

   /*
   * Destructor.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   FieldIoReal<D,RFRT,RFKT,FFTT>::~FieldIoReal()
   {}

   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   // %%%%%%% Dummy Functions : Not implemented in base class %%%%
   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   /*
   * These are virtual functions with a default implementation that 
   * should never be called, and which throws an Exception. An unusable 
   * implementation is provided rather than declaring these functions to
   * be pure virtual (= 0 suffix) because, in the absence of an pure
   * virtual functions, explicit instances of this template can be
   * compiled for the few use cases of interest (1, 2 or 3 dimensions,
   * with fields defined on the Cpu or Gpu). 
   */

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsRGrid(
                              std::istream &in,
                              DArray<RFRT >& fields,
                              UnitCell<D>& unitCell) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Read the data section of an array of fields in r-grid format.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsRGridData(
                              std::istream& in,
                              DArray< RFRT >& fields,
                              int nMonomer) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Read a single fields in r-grid format.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldRGrid(
                              std::istream &in,
                              RFRT & field,
                              UnitCell<D>& unitCell) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Write an array of fields in r-grid format
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsRGrid(
                              std::ostream &out,
                              DArray<RFRT > const & fields,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric,
                              bool writeMeshSize) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Read a single fields in r-grid format
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldRGrid(
                              std::ostream &out,
                              RFRT const & field,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Read an array of fields in k-grid format
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsKGrid(
                              std::istream &in,
                              DArray<RFKT >& fields,
                              UnitCell<D>& unitCell) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   template <int D, class RFRT, class RFKT, class FFTT>
   void
   FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsKGrid(
                              std::ostream &out,
                              DArray<RFKT > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFKT& out) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertKGridToBasis(
                              RFKT const & in,
                              DArray<double>& out,
                              bool checkSymmetry,
                              double epsilon) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Test if an real field DFT has the declared space group symmetry.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   bool FieldIoReal<D,RFRT,RFKT,FFTT>::hasSymmetry(
                              RFKT const & in, 
                              double epsilon,
                              bool verbose) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::replicateUnitCell(
                              std::ostream &out,
                              DArray< RFRT > const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   /*
   * Expand dimension of an array of r-grid fields, write to ostream.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::expandRGridDimension(
                              std::ostream &out,
                              DArray<RFRT> const & fields,
                              UnitCell<D> const & unitCell,
                              int d,
                              DArray<int> const& newGridDimensions) const
   {  UTIL_THROW("Unimplemented function in FieldIoReal base class"); }

   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   // %%%%%%% Shared Functions : Implemented in this template %%%%
   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   /*
   * Create associations with other members of parent Domain.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void
   FieldIoReal<D,RFRT,RFKT,FFTT>::associate(
                    Mesh<D> const & mesh,
                    FFTT const & fft,
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
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::setFileMaster(
                              FileMaster const & fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   // Field File IO - Symmetry-Adapted Basis Format

   /*
   * Read an array of fields in basis format from an input stream.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsBasis(
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

      // Check allocation of fields container, allocate if necessary
      if (fields.isAllocated()) {
         int nMonomerFields, fieldCapacity;
         inspectArrays(fields, fieldCapacity, nMonomerFields);
         UTIL_CHECK(nMonomerFields == nMonomer);
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(nBasisIn);
         }
      }

      // Read field data (components in basis)
      Prdc::readBasisData(in, fields, unitCell, mesh(), basis(), nBasisIn);
   }

   /*
   * Read a single field in basis format from an input stream
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldBasis(
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

      // Check that it only read 1 field
      UTIL_CHECK(fields.capacity() == 1);

      // Copy data from local array fields to function parameter field
      field = fields[0];
   }

   /*
   * Write an array of fields in basis format to an output stream.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsBasis(
                              std::ostream &out,
                              DArray< DArray<double> > const & fields,
                              UnitCell<D> const & unitCell) const
   {
      // Inspect fields to obtain nMonomer and fieldCapacity
      int nMonomer;
      int fieldCapacity;
      inspectArrays(fields, fieldCapacity, nMonomer);

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
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldBasis(
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

   // Field File IO - R-Grid Format
   // K-Grid Field Format

   // Field Format Conversion Functions - Symmetry Adapted Basis

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertBasisToKGrid(
                              DArray< DArray <double> > const & in,
                              DArray< RFKT >& out) const
   {
      // Inspect input and output field containers
      int nMonomer, nMonomerOut, capacity;
      inspectArrays(in, capacity, nMonomer);
      IntVec<D> dimensions;
      inspectFields(out, dimensions, nMonomerOut);
      UTIL_CHECK(nMonomer == nMonomerOut);

      // Convert fields for all monomer types
      for (int i = 0; i < nMonomer; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertKGridToBasis(
                              DArray< RFKT > const & in,
                              DArray< DArray <double> > & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      // Inspect input and output containers
      int nMonomer, nMonomerOut, capacity;
      IntVec<D> dimensions;
      inspectFields(in, dimensions, nMonomer);
      inspectArrays(out, capacity, nMonomerOut);
      UTIL_CHECK(nMonomer == nMonomerOut);

      // Convert fields for all monomer types
      for (int i = 0; i < nMonomer; ++i) {
         convertKGridToBasis(in[i], out[i], checkSymmetry, epsilon);
      }
   }

   /*
   * Test if an RFRT has declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   bool FieldIoReal<D,RFRT,RFKT,FFTT>::hasSymmetry(
                              RFRT const & in, 
                              double epsilon,
                              bool verbose) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_, epsilon, verbose);
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertBasisToRGrid(
                              DArray<double> const & in,
                              RFRT& out) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      convertBasisToKGrid(in, workDft_);
      fft().inverseTransformSafe(workDft_, out);
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertBasisToRGrid(
                              DArray< DArray <double> > const & in,
                              DArray< RFRT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkAllocateField(workDft_, mesh().dimensions());

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransformSafe(workDft_, out[i]);
      }
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertRGridToBasis(
                              RFRT const & in,
                              DArray<double> & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      convertKGridToBasis(workDft_, out, checkSymmetry, epsilon);
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertRGridToBasis(
                              DArray< RFRT > const & in,
                              DArray< DArray <double> > & out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      // Inspect input and output containers
      int nMonomer, nMonomerOut, capacity;
      IntVec<D> dimensions;
      inspectFields(in, dimensions, nMonomer);
      UTIL_CHECK(mesh().dimensions() == dimensions);
      inspectArrays(out, capacity, nMonomerOut);
      UTIL_CHECK(nMonomer == nMonomerOut);
      checkAllocateField(workDft_, dimensions);

      // Convert RGrid -> KGrid -> Basis for each field
      for (int i = 0; i < nMonomer; ++i) {
         fft().forwardTransform(in[i], workDft_);
         convertKGridToBasis(workDft_, out[i], checkSymmetry, epsilon);
      }
   }

   // FFT conversion functions (KGrid <-> RGrid) [Same in Rpc and Rpg]

   /*
   * Apply inverse FFT to an array of k-grid fields.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertKGridToRGrid(
                              DArray< RFKT > & in,
                              DArray< RFRT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   /*
   * Apply inverse FFT to a single k-grid field.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertKGridToRGrid(
                              RFKT& in, RFRT& out) const
   {
      fft().inverseTransformSafe(in, out);
   }

   /*
   * Apply forward FFT to an array of r-grid fields.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertRGridToKGrid(
                              DArray< RFRT > const & in,
                              DArray< RFKT >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Apply forward FFT to a single r-grid field.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::convertRGridToKGrid(
                              RFRT const & in,
                              RFKT& out) const
   {
      fft().forwardTransform(in, out);
   }

   // Grid Manipulation Utilities

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
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsBasis(
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
   * Open-close a file and read single fields in basis format.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldBasis(
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
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsBasis(
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
   * Write a single field in basis format, open and close the file.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldBasis( 
                              std::string filename,
                              DArray<double> const & field,
                              UnitCell<D> const & unitCell) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldBasis(file, field, unitCell);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsRGrid(
                              std::string filename,
                              DArray< RFRT >& fields,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldRGrid(
                              std::string filename,
                              RFRT & field,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldRGrid(file, field, unitCell);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsRGrid(
                              std::string filename,
                              DArray< RFRT > const & fields,
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

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldRGrid(
                              std::string filename,
                              RFRT const & field,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field, unitCell, isSymmetric);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldsKGrid(
                              std::string filename,
                              DArray< RFKT >& fields,
                              UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldsKGrid(
                              std::string filename,
                              DArray< RFKT > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell, isSymmetric);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::replicateUnitCell(
                              std::string filename,
                              DArray<RFRT> const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      replicateUnitCell(file, fields, unitCell, replicas);
      file.close();
   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::expandRGridDimension(
                              std::string filename,
                              DArray<RFRT> const & fields,
                              UnitCell<D> const & unitCell, int d,
                              DArray<int> newGridDimensions) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      expandRGridDimension(file, fields, unitCell, d, newGridDimensions);
      file.close();
   }

   // File Header IO Utilities

   /*
   * Read common part of field header and extract the number of monomers
   * (i.e., number of fields) and unitCell from the file.
   */
   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::readFieldHeader(
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

      // Check for presence of group name
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

      }

   }

   template <int D, class RFRT, class RFKT, class FFTT>
   void FieldIoReal<D,RFRT,RFKT,FFTT>::writeFieldHeader(
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

} // namespace Prdc
} // namespace Pscf
#endif
