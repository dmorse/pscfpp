#ifndef RPC_FIELD_IO_TPP
#define RPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <prdc/field/FieldIoReal.tpp>
#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>

#include <prdc/field/fieldIoUtil.h>
#include <prdc/crystal/fieldHeader.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>
#include <pscf/mesh/Mesh.h>

#if 0
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/SpaceGroup.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>

#include <util/misc/Log.h>
#endif

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

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
      checkAllocateFields(fields, mesh().dimensions(), nMonomer);

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, mesh().dimensions(), nMonomer);
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
      checkAllocateFields(fields, mesh().dimensions(), nMonomer);
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, mesh().dimensions(), nMonomer);
      // Rpg:: Copy host -> device
   }

   /*
   * Read a single fields in r-grid format.
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
   * Write an array of fields in r-grid format
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
      inspectFields(fields, meshDimensions, nMonomer);

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
      Prdc::writeRGridData(out, fields, meshDimensions, nMonomer);
   }

   /*
   * Read a single fields in r-grid format
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
      // Rpg:: Allocate host arrays
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

      checkAllocateFields(fields, mesh().dimensions(), nMonomer);

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readKGridData(in, fields, 
                          fields[0].dftDimensions(), nMonomer);
      // Rpg:: Copy host -> device
   }

   template <int D>
   void
   FieldIo<D>::writeFieldsKGrid(
                              std::ostream &out,
                              DArray<RFieldDft<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              bool isSymmetric) const
   {
      // Inspect fields array - infer nMonomer and dimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, meshDimensions, nMonomer);
      IntVec<D> dftDimensions = fields[0].dftDimensions();

      // Write file
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeMeshDimensions(out, meshDimensions);

      // Write data
      // Rpg:: Allocate host arrays
      // Rpg:: Copy device -> host
      Prdc::writeKGridData(out, fields, dftDimensions, nMonomer);
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(
                              DArray<double> const & in,
                              RFieldDft<D>& out) const
   {
      // Rpg: Allocate host arrays  
      Prdc::convertBasisToKGrid(in, out, basis(), out.dftDimensions()); 
      // Rpg: Copy host -> device
   }

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
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
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
      inspectFields(fields, meshDimensions, nMonomer);

      // Rpg: Allocate hostArrays
      // Rpg: Copy device -> host
      Prdc::replicateUnitCell(out, fields, meshDimensions,
                              unitCell, replicas);
   }

   /*
   * Expand dimension of an array of r-grid fields, write to ostream.
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
