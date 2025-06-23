#ifndef RPG_FIELD_IO_TPP
#define RPG_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <rpg/field/HostDArrayComplex.h>

#include <pscf/math/complex.h>

#include <prdc/field/FieldIoReal.tpp>      // base class implementation
#include <prdc/field/fieldIoUtil.h>
#include <prdc/field/fieldArrayUtil.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/cuda/RFieldComparison.h>
#include <prdc/cuda/types.h>
#include <prdc/cuda/complex.h>
#include <prdc/cuda/resources.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/HostDArray.h>


namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : FieldIoReal<D, RField<D>, RFieldDft<D>, FFT<D> >()
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldIo<D>::~FieldIo()
   {}

   // Field Io in r-grid format

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   bool FieldIo<D>::readFieldsRGrid(
                              std::istream &in,
                              DArray<RField<D> >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header and check fields dimensions
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Allocate host arrays
      DArray< HostDArray<cudaReal> > hostFields;
      allocateArrays(hostFields, nMonomer, mesh().size());

      // Read data
      Prdc::readRGridData(in, hostFields, nMonomer, mesh().dimensions());

      // Copy device <- host
      copyArrays(fields, hostFields);

      // Return true iff the header contains a space group declaration
      return isSymmetric;
   }

   /*
   * Read the data section of an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::readFieldsRGridData(
                              std::istream& in,
                              DArray< RField<D> >& fields,
                              int nMonomer) const
   {
      // Precondition: Check dimensions of fields
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Allocate host arrays
      DArray< HostDArray<cudaReal> > hostFields;
      allocateArrays(hostFields, nMonomer, mesh().size());

      // Read data section of file
      Prdc::readRGridData(in, hostFields, nMonomer, mesh().dimensions());

      // Copy device <- host
      copyArrays(fields, hostFields);
   }

   /*
   * Read a single field in r-grid format.
   */
   template <int D>
   bool FieldIo<D>::readFieldRGrid(
                              std::istream &in,
                              RField<D> & field,
                              UnitCell<D>& unitCell) const
   {

      // Read header and check field dimensions
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(nMonomer == 1);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateField(field, mesh().dimensions());

      // Allocate host field
      HostDArray<cudaReal> hostField;
      hostField.allocate(mesh().size());

      // Read data section with one field
      Prdc::readRGridData(in, hostField, mesh().dimensions());

      // Copy device <- host
      field = hostField;

      // Return true iff the header contains a space group declaration
      return isSymmetric;
   }

   /*
   * Write an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::writeFieldsRGrid(
                              std::ostream &out,
                              DArray<RField<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric,
                              bool writeMeshSize) const
   {
      // Inspect fields array, check field dimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      UTIL_CHECK(meshDimensions == mesh().dimensions());
      int meshSize = mesh().size();
      UTIL_CHECK(fields[0].capacity() == meshSize);

      // Write header
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      }
      if (writeMeshSize){
         writeMeshDimensions(out, meshDimensions);
      }

      // Copy field data to host container
      DArray< HostDArray<cudaReal> > hostFields;
      allocateArrays(hostFields, nMonomer, meshSize);
      copyArrays(hostFields, fields);

      // Write data section
      Prdc::writeRGridData(out, hostFields, nMonomer, meshDimensions);
   }

   /*
   * Write a single field in r-grid format.
   */
   template <int D>
   void FieldIo<D>::writeFieldRGrid(
                              std::ostream &out,
                              RField<D> const & field,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric) const
   {
      IntVec<D> meshDimensions = field.meshDimensions();
      int meshSize = field.capacity();
      UTIL_CHECK(meshDimensions == mesh().dimensions());
      UTIL_CHECK(meshSize == mesh().size());

      // Write header
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell, isSymmetric);
         writeMeshDimensions(out, meshDimensions);
      }

      // Copy field (device) to hostField
      HostDArray<cudaReal> hostField;
      hostField.allocate(meshSize);
      hostField = field;

      // Write data from hostField
      Prdc::writeRGridData(out, hostField, meshDimensions);
   }

   // Field IO in k-grid format

   /*
   * Read an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D>::readFieldsKGrid(
                           std::istream &in,
                           DArray<RFieldDft<D> >& fields,
                           UnitCell<D>& unitCell) const
   {
      // Read header and validate field mesh dimensions
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());
      IntVec<D> dftDimensions = fields[0].dftDimensions();
      int capacity = fields[0].capacity();

      // Allocate hostFields
      DArray< HostDArrayComplex > hostFields;
      allocateArrays(hostFields, nMonomer, capacity);

      // Read data into hostFields
      Prdc::readKGridData(in, hostFields, nMonomer, dftDimensions);

      // Copy device <- host
      copyArrays(fields, hostFields);
   }

   /*
   * Write an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D>::writeFieldsKGrid(
                              std::ostream &out,
                              DArray<RFieldDft<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      // Read header and validate field mesh dimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      UTIL_CHECK(mesh().dimensions() == meshDimensions);
      IntVec<D> dftDimensions = fields[0].dftDimensions();
      int capacity = fields[0].capacity();

      // Write header
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeMeshDimensions(out, meshDimensions);

      // Copy data from device to hostFields
      DArray< HostDArrayComplex > hostFields;
      allocateArrays(hostFields, nMonomer, capacity);
      copyArrays(hostFields, fields);

      // Write data from hostFields
      Prdc::writeKGridData(out, hostFields, nMonomer, dftDimensions);
   }

   /*
   * Convert a single field from basis to k-grid format.
   */
   template <int D>
   void FieldIo<D>::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D>& out) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() > 0);
      UTIL_CHECK(out.meshDimensions() == mesh().dimensions());

      // Allocate hostField
      HostDArrayComplex hostField;
      hostField.allocate(out.capacity());

      // Convert basis to k-grid on hostField
      Prdc::convertBasisToKGrid(in, hostField, basis(),
                                out.dftDimensions());

      // Copy out (device) <- host
      out = hostField;
   }

   /*
   * Write an array of fields from k-grid to basis format.
   */
   template <int D>
   void FieldIo<D>::convertKGridToBasis(
                              RFieldDft<D> const & in,
                              DArray<double>& out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.meshDimensions() == mesh().dimensions());
      UTIL_CHECK(out.capacity() > 0);

      // Copy k-grid input to hostField
      HostDArrayComplex hostField;
      hostField.allocate(in.capacity());
      hostField = in;

      // Convert k-grid host field to basis format
      Prdc::convertKGridToBasis(hostField, out, basis(),
                                in.dftDimensions(),
                                checkSymmetry, epsilon);
   }

   /*
   * Test if an real field DFT has the declared space group symmetry.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(
                              RFieldDft<D> const & in,
                              double epsilon,
                              bool verbose) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(in.meshDimensions() == mesh().dimensions());

      // Copy k-grid input to hostField
      HostDArrayComplex hostField;
      hostField.allocate(in.capacity());
      hostField = in;

      // Check symmetry of hostField
      return Prdc::hasSymmetry(hostField, basis(), in.dftDimensions(),
                               epsilon, verbose);
   }

   /*
   * Compare two fields in r-grid format, output report to Log file.
   */
   template <int D>
   void FieldIo<D>::compareFieldsRGrid(DArray< RField<D> > const & field1,
                                       DArray< RField<D> > const & field2) 
   const
   {
      RFieldComparison<D> comparison;
      comparison.compare(field1, field2);

      Log::file() << "\n Real-space field comparison results"
                  << std::endl;
      Log::file() << "     Maximum Absolute Difference:   "
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   "
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   /*
   * Multiply a field in r-grid format by a constant factor. 
   */
   template <int D>
   void FieldIo<D>::scaleFieldRGrid(
                              RField<D> & field,
                              double factor) const
   {
      UTIL_CHECK(field.isAllocated());
      VecOp::mulEqS(field, factor);
   }

   /*
   * Replicate the unit cell for an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D>::replicateUnitCell(
                              std::ostream &out,
                              DArray< RField<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const

   {
      // Inspect fields to obtain nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      UTIL_CHECK(meshDimensions == mesh().dimensions());
      int capacity = fields[0].capacity();

      // Copy k-grid input to hostField
      DArray< HostDArray<cudaReal> > hostFields;
      allocateArrays(hostFields, nMonomer, capacity);
      copyArrays(hostFields, fields);

      Prdc::replicateUnitCell(out, hostFields, meshDimensions,
                              unitCell, replicas);
   }

   /*
   * Expand spatial dimension of an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D>::expandRGridDimension(
                              std::ostream &out,
                              DArray< RField<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              int d,
                              DArray<int> const& newGridDimensions) const
   {
      // Inspect fields to obtain nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      UTIL_CHECK(meshDimensions == mesh().dimensions());
      int capacity = fields[0].capacity();

      // Copy k-grid input fields to hostFields
      DArray< HostDArray<cudaReal> > hostFields;
      allocateArrays(hostFields, nMonomer, capacity);
      copyArrays(hostFields, fields);

      Prdc::expandRGridDimension(out, hostFields, meshDimensions,
                                 unitCell, d, newGridDimensions);
   }

}
}
#endif
