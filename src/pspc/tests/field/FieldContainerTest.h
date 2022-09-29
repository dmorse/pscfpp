#ifndef PSPC_FIELD_CONTAINER_TEST_H
#define PSPC_FIELD_CONTAINER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/RFieldComparison.h>

#include <pspc/field/FieldContainer.h>
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

class FieldContainerTest : public UnitTest 
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
                       DArray< DArray<double> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(nBasis);
      }
   }

   // Allocate an array of r-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
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

   template <int D>
   void writeFields(std::string filename, Domain<D>& domain,
                   DArray< RField<D> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsRGrid(out, fields, domain.unitCell());
      out.close();
   }

   template <int D>
   void writeFields(std::string filename, Domain<D>& domain,
                   DArray< RFieldDft<D> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsKGrid(out, fields, domain.unitCell());
      out.close();
   }

   void testSetBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      //std::cout << "\n domain.basis().nBasis() = "
      //          << domain.basis().nBasis() << "\n";

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      FieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocated());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      fields.setBasis(bf);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      //std::cout << comparison.maxDiff() << std::endl;
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      domain.fieldIo().convertRGridToBasis(fields.rgrid(), bf_1);
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testSetRGrid_1_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      readFields("in/w_bcc.rf", domain, rf);

      FieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocated());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      fields.setRGrid(rf);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());

      RFieldComparison<3> comparison;
      comparison.compare(rf, fields.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testSetRGrid_2_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      DArray< RField<3> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      domain.fieldIo().convertBasisToRGrid(bf, rf);

      FieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocated());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      bool isSymmetric = true;
      fields.setRGrid(rf, isSymmetric);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }
};

TEST_BEGIN(FieldContainerTest)
TEST_ADD(FieldContainerTest, testSetBasis_bcc)
TEST_ADD(FieldContainerTest, testSetRGrid_1_bcc)
TEST_ADD(FieldContainerTest, testSetRGrid_2_bcc)
TEST_END(FieldContainerTest)

#endif
