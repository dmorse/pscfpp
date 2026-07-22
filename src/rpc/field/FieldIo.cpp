/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"                 // class header

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldComparison.h>
#include <pscf/cpu/complex.h>

#include <rp/field/FieldIoBase.tpp>  // base class template implementation

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   bool FieldIo<D, CppTp<D> >::readFieldsRGrid(
                              std::istream &in,
                              DArray< RField<D, CppTp<D> > >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      bool isSymmetric;
      FieldIoBase::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, FieldIoBase::mesh().dimensions());
      checkAllocateFields(fields, nMonomer, FieldIoBase::mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, nMonomer, FieldIoBase::mesh().dimensions());
      // Rpg:: Copy host -> device

      return isSymmetric;
   }

   /*
   * Read the data section of an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::readFieldsRGridData(
                              std::istream& in,
                              DArray< RField<D, CppTp<D> > >& fields,
                              int nMonomer) const
   {
      checkAllocateFields(fields, nMonomer, FieldIoBase::mesh().dimensions());
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, nMonomer, 
                          FieldIoBase::mesh().dimensions());
      // Rpg:: Copy host -> device
   }

   /*
   * Read a single field in r-grid format.
   */
   template <int D>
   bool FieldIo<D, CppTp<D> >::readFieldRGrid(
                              std::istream &in,
                              RField<D, CppTp<D> > & field,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      bool isSymmetric;
      FieldIoBase::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(nMonomer == 1);
      readMeshDimensions(in, FieldIoBase::mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      checkAllocateField(field, FieldIoBase::mesh().dimensions());
      Prdc::readRGridData(in, field, FieldIoBase::mesh().dimensions());
      // Rpg:: Copy host -> device

      return isSymmetric;
   }

   /*
   * Write an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::writeFieldsRGrid(
                              std::ostream &out,
                              DArray< RField<D, CppTp<D> > > const & fields,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric,
                              bool writeMeshSize) const
   {
      // Inspect fields array, infer nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);

      // Write header
      if (writeHeader){
         FieldIoBase::writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      }
      if (writeMeshSize){
         writeMeshDimensions(out, meshDimensions);
      }

      // Write data section
      // Rpg:: Allocate host arrays
      // Rpg:: Copy device -> host
      Prdc::writeRGridData(out, fields, nMonomer, meshDimensions);
   }

   /*
   * Write a single field in r-grid format.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::writeFieldRGrid(
                              std::ostream &out,
                              RField<D, CppTp<D> > const & field,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric) const
   {
      IntVec<D> meshDimensions = field.meshDimensions();

      // Write header
      if (writeHeader) {
         FieldIoBase::writeFieldHeader(out, 1, unitCell, isSymmetric);
         writeMeshDimensions(out, meshDimensions);
      }

      // Write data
      // Rpg:: Allocate host array
      // Rpg:: Copy device -> host
      Prdc::writeRGridData(out, field, meshDimensions);
   }

   /*
   * Read an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::readFieldsKGrid(
                              std::istream &in,
                              DArray< RFieldDft<D, CppTp<D> > >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      bool isSymmetric;
      FieldIoBase::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, FieldIoBase::mesh().dimensions());

      checkAllocateFields(fields, nMonomer, FieldIoBase::mesh().dimensions());
      IntVec<D> dftDimensions = fields[0].dftDimensions();

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readKGridData(in, fields, nMonomer, dftDimensions);
      // Rpg:: Copy host -> device
   }

   /*
   * Write an array of fields in k-grid format
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::writeFieldsKGrid(
                              std::ostream &out,
                              DArray< RFieldDft<D, CppTp<D> > > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      // Inspect fields array - infer nMonomer and dimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      IntVec<D> dftDimensions = fields[0].dftDimensions();

      // Write file
      FieldIoBase::writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeMeshDimensions(out, meshDimensions);

      // Write data
      // Rpg:: Allocate host arrays
      // Rpg:: Copy device -> host
      Prdc::writeKGridData(out, fields, nMonomer, dftDimensions);
   }

   /*
   * Convert an array of fields from basis to k-grid format.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D, CppTp<D> >& out) const
   {
      // Rpg: Allocate host array
      Prdc::convertBasisToKGrid(in, out, FieldIoBase::basis(), out.dftDimensions());
      // Rpg: Copy host -> device
   }

   /*
   * Convert an array of fields from k-grid to basis format.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::convertKGridToBasis(
                          RFieldDft<D, CppTp<D> > const & in,
                          DArray<double>& out,
                          bool checkSymmetry,
                          double epsilon) const
   {
      // Rpg: Allocate host arrays
      // Rpg: Copy device -> host
      Prdc::convertKGridToBasis(in, out, FieldIoBase::basis(), 
                          in.dftDimensions(),
                          checkSymmetry, epsilon);
   }

   /*
   * Test if an real field DFT has the declared space group symmetry.
   */
   template <int D>
   bool FieldIo<D, CppTp<D> >::hasSymmetry(
                          RFieldDft<D, CppTp<D> > const & in,
                          double epsilon,
                          bool verbose) const
   {
      // Rpg:: Allocate host array
      // Rpg: Copy device -> host
      return Prdc::hasSymmetry(in, FieldIoBase::basis(), in.dftDimensions(),
                               epsilon, verbose);
   }

   /*
   * Compare two fields in r-grid format, output report to Log file.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::compareFieldsRGrid(
                          DArray< RField<D, CppTp<D> > > const & field1,
                          DArray< RField<D, CppTp<D> > > const & field2) const
   {
      RFieldComparison<D, CppTp<D> > comparison;
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
   void FieldIo<D, CppTp<D> >::scaleFieldRGrid(
                              RField<D, CppTp<D> > & field,
                              double factor) const
   {
      int n = field.capacity();
      for (int i = 0; i < n; ++i) {
         field[i] *= factor;
      }
   }

   /*
   * Replicate the unit cell for an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::replicateUnitCell(
                              std::ostream &out,
                              DArray< RField<D, CppTp<D> > > const & fields,
                              UnitCell<D> const & unitCell,
                              IntVec<D> const & replicas) const

   {
      // Inspect fields to obtain nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);

      // Rpg: Allocate hostArrays
      // Rpg: Copy device -> host
      Prdc::replicateUnitCell(out, fields, meshDimensions,
                              unitCell, replicas);
   }

   /*
   * Expand spatial dimension of an array of r-grid fields.
   */
   template <int D>
   void FieldIo<D, CppTp<D> >::expandRGridDimension(
                              std::ostream &out,
                              DArray< RField<D, CppTp<D> > > const & fields,
                              UnitCell<D> const & unitCell,
                              int d,
                              DArray<int> const& newGridDimensions) const
   {
      IntVec<D> meshDimensions = fields[0].meshDimensions();

      // Rpg: Allocate hostArrays
      // Rpg: Copy device -> host
      Prdc::expandRGridDimension(out, fields, meshDimensions,
                                 unitCell, d, newGridDimensions);
   }

} // namespace Rp
} // namespace Pscf

// Explicit specialization definitions
namespace Pscf {
   namespace Rp {
      template class Rp::FieldIoBase<1, CppTp<1> >;
      template class Rp::FieldIoBase<2, CppTp<2> >;
      template class Rp::FieldIoBase<3, CppTp<3> >;
      template class Rp::FieldIo<1, CppTp<1> >;
      template class Rp::FieldIo<2, CppTp<2> >;
      template class Rp::FieldIo<3, CppTp<3> >;
   }
}
