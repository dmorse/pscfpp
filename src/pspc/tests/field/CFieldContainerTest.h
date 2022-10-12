#ifndef PSPC_C_FIELD_CONTAINER_TEST_H
#define PSPC_C_FIELD_CONTAINER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/RFieldComparison.h>

#include <pspc/field/CFieldContainer.h>
#include <pspc/field/Domain.h>
#include <pspc/field/FieldIo.h>
#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>

#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class CFieldContainerTest : public UnitTest 
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
      domain.readFieldHeader(in, nMonomer_);
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

      CFieldContainer<3> fields;
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
      std::cout << "\n nMonomer_      = " << nMonomer_ << "\n";
      std::cout << "\n rgrid capacity = " << fields.rgrid().capacity() << "\n";
      std::cout << "\n basis capacity = " << fields.rgrid().capacity() << "\n";
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
   }

};

TEST_BEGIN(CFieldContainerTest)
TEST_ADD(CFieldContainerTest, testAllocate)
TEST_END(CFieldContainerTest)

#endif
