#ifndef RPC_FIELD_IO_TPP
#define RPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"
#include <prdc/field/FieldIoReal.tpp>     // base class implementation 

#include <rpc/field/Domain.h>

#include <prdc/field/fieldIoUtil.h>
#include <prdc/crystal/fieldHeader.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/cpu/complex.h>
#include <prdc/cpu/types.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : FieldIoReal< D, RField<D>, RFieldDft<D>, FFT<D> >()
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldIo<D>::~FieldIo()
   {}

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::readFieldsRGrid(
                              std::istream &in,
                              DArray<RField<D> >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, nMonomer, mesh().dimensions());
      // Rpg:: Copy host -> device

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
      checkAllocateFields(fields, nMonomer, mesh().dimensions());
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, nMonomer, mesh().dimensions());
      // Rpg:: Copy host -> device
   }

   /*
   * Read a single field in r-grid format.
   */
   template <int D>
   void FieldIo<D>::readFieldRGrid(
                              std::istream &in,
                              RField<D> & field,
                              UnitCell<D>& unitCell) const
   {

      // Read header
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(nMonomer == 1);
      readMeshDimensions(in, mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      checkAllocateField(field, mesh().dimensions());
      Prdc::readRGridData(in, field, mesh().dimensions());
      // Rpg:: Copy host -> device
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
      // Inspect fields array, infer nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);

      // Write header
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
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
   void FieldIo<D>::writeFieldRGrid(
                              std::ostream &out,
                              RField<D> const & field,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool isSymmetric) const
   {
      IntVec<D> meshDimensions = field.meshDimensions();

      // Write header
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell, isSymmetric);
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
   void FieldIo<D>::readFieldsKGrid(
                              std::istream &in,
                              DArray<RFieldDft<D> >& fields,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      bool isSymmetric;
      readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in, mesh().dimensions());

      checkAllocateFields(fields, nMonomer, mesh().dimensions());
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
   void FieldIo<D>::writeFieldsKGrid(
                              std::ostream &out,
                              DArray<RFieldDft<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      // Inspect fields array - infer nMonomer and dimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);
      IntVec<D> dftDimensions = fields[0].dftDimensions();

      // Write file
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
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
   void FieldIo<D>::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D>& out) const
   {
      // Rpg: Allocate host array 
      Prdc::convertBasisToKGrid(in, out, basis(), out.dftDimensions()); 
      // Rpg: Copy host -> device
   }

   /*
   * Convert an array of fields from k-grid to basis format.
   */
   template <int D>
   void FieldIo<D>::convertKGridToBasis(
                              RFieldDft<D> const & in,
                              DArray<double>& out,
                              bool checkSymmetry,
                              double epsilon) const
   {
      // Rpg: Allocate host arrays
      // Rpg: Copy device -> host 
      Prdc::convertKGridToBasis(in, out, basis(), in.dftDimensions(),
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
      // Rpg:: Allocate host array
      // Rpg: Copy device -> host 
      return Prdc::hasSymmetry(in, basis(), in.dftDimensions(),
                               epsilon, verbose);
   }

   /*
   * Test if an real field DFT has the declared space group symmetry.
   */
   template <int D>
   void FieldIo<D>::scaleFieldBasis(
                              DArray<double> & field, 
                              double factor) const
   {
      int n = field.capacity();
      for (int i = 0; i < n; ++i) {
         field[i] *= factor;
      }
   }

   /*
   * Test if an real field DFT has the declared space group symmetry.
   */
   template <int D>
   void FieldIo<D>::scaleFieldRGrid(
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

      // Rpg: Allocate hostArrays
      // Rpg: Copy device -> host
      Prdc::replicateUnitCell(out, fields, meshDimensions,
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
      IntVec<D> meshDimensions = fields[0].meshDimensions();

      // Rpg: Allocate hostArrays
      // Rpg: Copy device -> host
      Prdc::expandRGridDimension(out, fields, meshDimensions,
                                 unitCell, d, newGridDimensions);
   }

}
}
#endif
