/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo_cp.h"              // partial specialization header

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldComparison.h>
#include <pscf/cpu/complex.h>

#include <rp/field/FieldIoBase.tpp>  // base class implementation

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   bool FieldIo<D, Rpc::Types<D> >::readFieldsRGrid(
                              std::istream &in,
                              DArray< RField<D> >& fields,
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
   void FieldIo<D, Rpc::Types<D> >::readFieldsRGridData(
                              std::istream& in,
                              DArray< RField<D> >& fields,
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
   bool FieldIo<D, Rpc::Types<D> >::readFieldRGrid(
                              std::istream &in,
                              RField<D> & field,
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
   void FieldIo<D, Rpc::Types<D> >::writeFieldsRGrid(
                              std::ostream &out,
                              DArray< RField<D> > const & fields,
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
   void FieldIo<D, Rpc::Types<D> >::writeFieldRGrid(
                              std::ostream &out,
                              RField<D> const & field,
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
   void FieldIo<D, Rpc::Types<D> >::readFieldsKGrid(
                              std::istream &in,
                              DArray< RFieldDft<D> >& fields,
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
   void FieldIo<D, Rpc::Types<D> >::writeFieldsKGrid(
                              std::ostream &out,
                              DArray< RFieldDft<D> > const & fields,
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
   void FieldIo<D, Rpc::Types<D> >::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D>& out) const
   {
      // Rpg: Allocate host array
      Prdc::convertBasisToKGrid(in, out, FieldIoBase::basis(), out.dftDimensions());
      // Rpg: Copy host -> device
   }

   /*
   * Convert an array of fields from k-grid to basis format.
   */
   template <int D>
   void FieldIo<D, Rpc::Types<D> >::convertKGridToBasis(
                          RFieldDft<D> const & in,
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
   bool FieldIo<D, Rpc::Types<D> >::hasSymmetry(
                          RFieldDft<D> const & in,
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
   void FieldIo<D, Rpc::Types<D> >::compareFieldsRGrid(
                          DArray< RField<D> > const & field1,
                          DArray< RField<D> > const & field2) const
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
   void FieldIo<D, Rpc::Types<D> >::scaleFieldRGrid(
                              RField<D> & field,
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
   void FieldIo<D, Rpc::Types<D> >::replicateUnitCell(
                              std::ostream &out,
                              DArray< RField<D> > const & fields,
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
   void FieldIo<D, Rpc::Types<D> >::expandRGridDimension(
                              std::ostream &out,
                              DArray< RField<D> > const & fields,
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
      template class Rp::FieldIoBase<1, Rpc::Types<1> >;
      template class Rp::FieldIoBase<2, Rpc::Types<2> >;
      template class Rp::FieldIoBase<3, Rpc::Types<3> >;
      template class Rp::FieldIo<1, Rpc::Types<1> >;
      template class Rp::FieldIo<2, Rpc::Types<2> >;
      template class Rp::FieldIo<3, Rpc::Types<3> >;
   }
}
