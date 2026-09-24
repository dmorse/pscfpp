/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/fieldIo/fieldCheck.h>
#include <pscf/backend/cuda/ConstHostArray.h>
#include <pscf/backend/cuda/HostArray.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/complex.h>

#include <rp/field/FieldIoBase.tpp>   // base class implementation
#include <rp/field/FieldIo_u.h>       // class header

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   // Field Io in r-grid format

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   bool FieldIo<D,CUT>::readFieldsRGrid(
                              std::istream &in,
                              DArray< RField<D,CUT> >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header and check fields dimensions
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Setup host arrays
      DArray< HostArray<RealT,CUT> > hostFields;
      associateArrays(hostFields, fields);

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
   void FieldIo<D,CUT>::readFieldsRGridData(
                              std::istream& in,
                              DArray< RField<D,CUT> >& fields,
                              int nMonomer) const
   {
      // Precondition: Check dimensions of fields
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Allocate host arrays
      DArray< HostArray<RealT,CUT> > hostFields;
      associateArrays(hostFields, fields);

      // Read data section of file
      Prdc::readRGridData(in, hostFields, nMonomer, mesh().dimensions());

      // Copy device <- host
      copyArrays(fields, hostFields);
   }

   /*
   * Read a single field in r-grid format.
   */
   template <int D>
   bool FieldIo<D,CUT>::readFieldRGrid(
                              std::istream &in,
                              RField<D,CUT> & field,
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
      HostArray<RealT,CUT> hostField;
      hostField.associate(field);

      // Read data section with one field
      Prdc::readRGridData(in, hostField, mesh().dimensions());

      // Copy to device from host
      field = hostField;

      // Return true iff the header contains a space group declaration
      return isSymmetric;
   }

   /*
   * Write an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D,CUT>::writeFieldsRGrid(
                              std::ostream &out,
                              DArray< RField<D,CUT> > const & fields,
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
      DArray< ConstHostArray<RealT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Write data section
      Prdc::writeRGridData(out, hostFields, nMonomer, meshDimensions);
   }

   /*
   * Write a single field in r-grid format.
   */
   template <int D>
   void FieldIo<D,CUT>::writeFieldRGrid(
                              std::ostream &out,
                              RField<D,CUT> const & field,
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
      ConstHostArray<RealT,CUT> hostField;
      hostField = field;

      // Write data from hostField
      Prdc::writeRGridData(out, hostField, meshDimensions);
   }

   // Field IO in k-grid format

   /*
   * Read an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D,CUT>::readFieldsKGrid(
                           std::istream &in,
                           DArray< RFieldDft<D,CUT> >& fields,
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
      DArray< HostArray<ComplexT,CUT> > hostFields;
      associateArrays(hostFields, fields);

      // Read data into hostFields
      Prdc::readKGridData(in, hostFields, nMonomer, dftDimensions);

      // Copy device <- host
      copyArrays(fields, hostFields);
   }

   /*
   * Write an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D,CUT>::writeFieldsKGrid(
                              std::ostream &out,
                              DArray< RFieldDft<D,CUT> > const & fields,
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
      DArray< ConstHostArray<ComplexT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Write data from hostFields
      Prdc::writeKGridData(out, hostFields, nMonomer, dftDimensions);
   }

   /*
   * Convert a single field from basis to k-grid format.
   */
   template <int D>
   void FieldIo<D,CUT>::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D,CUT>& out) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() > 0);
      UTIL_CHECK(out.meshDimensions() == mesh().dimensions());

      // Allocate hostField
      HostArray<ComplexT,CUT> hostField;
      hostField.associate(out);

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
   void FieldIo<D,CUT>::convertKGridToBasis(
                              RFieldDft<D,CUT> const & in,
                              DArray<double>& out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.meshDimensions() == mesh().dimensions());
      UTIL_CHECK(out.capacity() > 0);

      // Copy k-grid input to hostField
      ConstHostArray<ComplexT,CUT> hostField;
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
   bool FieldIo<D,CUT>::hasSymmetry(
                              RFieldDft<D,CUT> const & in,
                              double epsilon,
                              bool verbose) const
   {
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(in.meshDimensions() == mesh().dimensions());

      // Copy k-grid input to hostField
      ConstHostArray<ComplexT,CUT> hostField;
      hostField = in;

      // Check symmetry of hostField
      return Prdc::hasSymmetry(hostField, basis(), in.dftDimensions(),
                               epsilon, verbose);
   }

   /*
   * Compare two fields in r-grid format, output report to Log file.
   */
   template <int D>
   void FieldIo<D,CUT>::compareFieldsRGrid(
                             DArray< RField<D,CUT> > const & field1,
                             DArray< RField<D,CUT> > const & field2) const
   {
      RFieldComparison<D,CUT> comparison;
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
   void FieldIo<D,CUT>::scaleFieldRGrid(
                              RField<D,CUT> & field,
                              double factor) const
   {
      UTIL_CHECK(field.isAllocated());
      VecOp::mulEqS(field, factor);
   }

   /*
   * Replicate the unit cell for an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D,CUT>::replicateUnitCell(
                              std::ostream &out,
                              DArray< RField<D,CUT> > const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const

   {
      // Inspect fields to obtain nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      UTIL_CHECK(meshDimensions == mesh().dimensions());
      int capacity = fields[0].capacity();

      // Copy r-grid input from device to host
      DArray< ConstHostArray<RealT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Compute replicated fields and write to a file
      Prdc::replicateUnitCell(out, hostFields, meshDimensions,
                              unitCell, replicas);
   }

   /*
   * Expand spatial dimension of an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D,CUT>::expandRGridDimension(
                              std::ostream &out,
                              DArray< RField<D,CUT> > const & fields,
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
      DArray< ConstHostArray<RealT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      Prdc::expandRGridDimension(out, hostFields, meshDimensions,
                                 unitCell, d, newGridDimensions);
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FieldIoBase<1,CUT>;
      template class FieldIoBase<2,CUT>;
      template class FieldIoBase<3,CUT>;
      template class FieldIo<1,CUT>;
      template class FieldIo<2,CUT>;
      template class FieldIo<3,CUT>;
   }
}
