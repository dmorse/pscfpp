#ifndef CPC_FIELD_IO_TPP
#define CPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Contains declarations needed to parse other templates
#include <pscf/cpu/complex.h>

#include "FieldIo.h"
#include <prdc/fieldIo/fieldCheck.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/field/cpu/CFieldComparison.h>
#include <prdc/fieldIo/cFieldIo.h>
#include <prdc/fieldIo/rFieldIo.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <cp/field/FieldIo.tpp>

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /*
   * Read an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::readFields(std::istream& in,
                               DArray< CField<D> >& fields,
                               UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      readFieldHeader(in, nMonomer, unitCell);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readCFieldsData(in, fields, mesh().dimensions());
      // Rpg:: Copy host -> device

   }

   /*
   * Read the data section of an array of complex fields.
   */
   template <int D>
   void FieldIo<D>::readFieldsData(std::istream& in,
                                   DArray< CField<D> >& fields,
                                   int nMonomer) const
   {
      checkAllocateFields(fields, nMonomer, mesh().dimensions());
      // Rpg:: Allocate host arrays
      Prdc::readCFieldsData(in, fields, mesh().dimensions());
      // Rpg:: Copy host -> device
   }

   /*
   * Read a single complex field.
   */
   template <int D>
   void FieldIo<D>::readField(std::istream &in,
                              CField<D> & field,
                              UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      readFieldHeader(in, nMonomer, unitCell);
      UTIL_CHECK(nMonomer == 1);
      readMeshDimensions(in, mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      checkAllocateField(field, mesh().dimensions());
      Prdc::readCFieldData(in, field, mesh().dimensions());
      // Rpg:: Copy host -> device

   }

   /*
   * Read an array of real-valued fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::istream& in,
                                    DArray< CField<D> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      // Read header
      int nMonomer;
      readFieldHeader(in, nMonomer, unitCell);
      readMeshDimensions(in, mesh().dimensions());
      checkAllocateFields(fields, nMonomer, mesh().dimensions());

      // Read data
      // Rpg:: Allocate host arrays
      Prdc::readRGridData(in, fields, nMonomer, mesh().dimensions());
      // Rpg:: Copy host -> device

   }

   /*
   * Write an array of fields in r-grid format.
   */
   template <int D>
   void FieldIo<D>::writeFields(
                              std::ostream &out,
                              DArray<CField<D> > const & fields,
                              UnitCell<D> const & unitCell,
                              bool writeHeader,
                              bool writeMeshSize) const
   {
      // Inspect fields array, infer nMonomer and meshDimensions
      int nMonomer;
      IntVec<D> meshDimensions;
      inspectFields(fields, nMonomer, meshDimensions);

      // Write header
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell);
      }
      if (writeMeshSize){
         writeMeshDimensions(out, meshDimensions);
      }

      // Write data section
      // Rpg:: Allocate host arrays
      // Rpg:: Copy device -> host
      using AT = CField<D>;
      using CT = typename CField<D>::ValueType;
      using RT = typename CField<D>::RealType;
      Prdc::writeCFieldsData<D,AT,CT,RT>(out, fields, meshDimensions);
   }

   /*
   * Write a single field in r-grid format.
   */
   template <int D>
   void FieldIo<D>::writeField(std::ostream &out,
                               CField<D> const & field,
                               UnitCell<D> const & unitCell,
                               bool writeHeader) const
   {
      IntVec<D> meshDimensions = field.meshDimensions();

      // Write header
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell);
         writeMeshDimensions(out, meshDimensions);
      }

      // Write data
      // Rpg:: Allocate host array
      // Rpg:: Copy device -> host
      using AT = CField<D>;
      using CT = typename CField<D>::ValueType;
      using RT = typename CField<D>::RealType;
      Prdc::writeCFieldData<D,AT,CT,RT>(out, field, meshDimensions);
   }

   #if 0
   /*
   * Compare two fields in r-grid format, output report to Log file.
   */
   template <int D>
   void FieldIo<D>::compareFields(DArray< CField<D> > const & field1,
                                  DArray< CField<D> > const & field2)
   const
   {
      CFieldComparison<D> comparison;
      comparison.compare(field1, field2);

      Log::file() << "\n Real-space field comparison results"
                  << std::endl;
      Log::file() << "     Maximum Absolute Difference:   "
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   "
                  << comparison.rmsDiff() << "\n" << std::endl;
   }
   #endif

}
}
#endif
