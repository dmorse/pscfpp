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

#if 0
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

      // Setup local host arrays
      DArray< HostArray<RealT,CUT> > hostFields;
      associateArrays(hostFields, fields);

      // Read data
      Prdc::readRGridData(in, hostFields, nMonomer, mesh().dimensions());

      // Copy host -> device 
      copyArrays(fields, hostFields);
      dissociateArrays(hostFields);

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

      // Setup local host arrays
      DArray< HostArray<RealT,CUT> > hostFields;
      associateArrays(hostFields, fields);

      // Read data section of file
      Prdc::readRGridData(in, hostFields, nMonomer, mesh().dimensions());

      // Copy host -> device 
      copyArrays(fields, hostFields);
      dissociateArrays(hostFields);
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

      // Setup local host array
      HostArray<RealT,CUT> hostField;
      hostField.associate(field);

      // Read data section with one field
      Prdc::readRGridData(in, hostField, mesh().dimensions());

      // Copy from host to device 
      field = hostField;
      hostField.dissociate();

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

      dissociateArrays(hostFields);
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

      // Copy field data to host container
      ConstHostArray<RealT,CUT> hostField;
      hostField = field;

      // Write data from hostField
      Prdc::writeRGridData(out, hostField, meshDimensions);

      hostField.dissociate();
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

      // Allocate hostFields
      DArray< HostArray<ComplexT,CUT> > hostFields;
      associateArrays(hostFields, fields);

      // Read data into hostFields
      Prdc::readKGridData(in, hostFields, nMonomer, dftDimensions);

      // Copy device <- host
      copyArrays(fields, hostFields);

      dissociateArrays(hostFields);
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

      // Write header
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeMeshDimensions(out, meshDimensions);

      // Copy data from device to host container
      DArray< ConstHostArray<ComplexT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Write data from host container
      Prdc::writeKGridData(out, hostFields, nMonomer, dftDimensions);

      dissociateArrays(hostFields);
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

      // Setup host container for k-grid data
      HostArray<ComplexT,CUT> hostField;
      hostField.associate(out);

      // Convert basis to k-grid on host
      Prdc::convertBasisToKGrid(in, hostField, basis(),
                                out.dftDimensions());

      // Copy from host to device
      out = hostField;

      hostField.dissociate();
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

      // Copy k-grid data from device to const host container
      ConstHostArray<ComplexT,CUT> hostField;
      hostField = in;

      // Convert from k-grid to basis format on host
      Prdc::convertKGridToBasis(hostField, out, basis(),
                                in.dftDimensions(),
                                checkSymmetry, epsilon);

      hostField.dissociate();
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

      // Copy k-grid data from device to const host container
      ConstHostArray<ComplexT,CUT> hostField;
      hostField = in;

      // Check symmetry of k-grid data on host, return result
      return Prdc::hasSymmetry(hostField, basis(), in.dftDimensions(),
                               epsilon, verbose);

      hostField.dissociate();
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

      // Copy r-grid input from device to host
      DArray< ConstHostArray<RealT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Compute replicated fields and write to a file
      Prdc::replicateUnitCell(out, hostFields, meshDimensions,
                              unitCell, replicas);

      dissociateArrays(hostFields);
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

      // Copy k-grid data from device to const host container
      DArray< ConstHostArray<RealT,CUT> > hostFields;
      copyArrays(hostFields, fields);

      // Write fields on higher-dimensional grid to file
      Prdc::expandRGridDimension(out, hostFields, meshDimensions,
                                 unitCell, d, newGridDimensions);

      dissociateArrays(hostFields);
   }

}
}
#endif

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
