#ifndef RPC_C_FIELD_CONTAINER_TEST_H
#define RPC_C_FIELD_CONTAINER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/field/CFields.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Rpc;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;

class CFieldsTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
   }

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   // Open and read parameter header to initialize Domain<D> system.
   template <int D>
   void readParam(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readParam(in);
      in.close();
   }

   // Open and read file header to initialize Domain<D> system.
   template <int D>
   void readHeader(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readRGridFieldHeader(in, nMonomer_);
      in.close();
   }

   // Allocate an array of fields in symmetry adapated format.
   void allocateFields(int nMonomer, int nBasis,
                       DArray< DArray<double> > & fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(nBasis);
      }
   }

   // Allocate an array of r-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> const & dimensions,
                       DArray< RField<D> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(dimensions);
      }
   }

   template <int D>
   void readFields(std::string filename, Domain<D>& domain,
                   DArray< DArray<double> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsBasis(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void readFields(std::string filename, Domain<D>& domain,
                   DArray< RField<D> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsRGrid(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void writeFields(std::string filename, Domain<D>& domain,
                   DArray< DArray<double> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsBasis(out, fields, domain.unitCell());
      out.close();
   }

   void testAllocate()
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      CFields<3> fields;
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
   }

};

TEST_BEGIN(CFieldsTest)
TEST_ADD(CFieldsTest, testAllocate)
TEST_END(CFieldsTest)

#endif
