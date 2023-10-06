#ifndef PSPC_W_FIELD_CONTAINER_TEST_H
#define PSPC_W_FIELD_CONTAINER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>


#include <pspc/field/WFieldContainer.h>
#include <pspc/field/Domain.h>
#include <pspc/field/FieldIo.h>

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
using namespace Pscf::Pspc;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;

class WFieldContainerTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      openLogFile("out/wFieldContainerTestLogFile");
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

   void testAllocate_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);
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

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
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

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
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

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedRGrid());
      TEST_ASSERT(fields.isAllocatedBasis());
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

   void testReadBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedBasis());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      std::ifstream in;
      openInputFile("in/w_bcc.bf", in);
      fields.readBasis(in, domain.unitCell());
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      domain.fieldIo().convertRGridToBasis(fields.rgrid(), bf_1);
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testReadRGrid_1_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      readFields("in/w_bcc.rf", domain, rf);

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedBasis());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      std::ifstream in;
      openInputFile("in/w_bcc.rf", in);
      fields.readRGrid(in, domain.unitCell());
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());

      RFieldComparison<3> comparison;
      comparison.compare(rf, fields.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testReadRGrid_2_bcc() 
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

      WFieldContainer<3> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.allocate(nMonomer_, domain.basis().nBasis(),
                      domain.mesh().dimensions());
      TEST_ASSERT(fields.isAllocatedBasis());
      TEST_ASSERT(!fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());
      TEST_ASSERT(fields.basis().capacity() == nMonomer_);
      TEST_ASSERT(fields.rgrid().capacity() == nMonomer_);

      std::ifstream in;
      openInputFile("in/w_bcc.rf", in);
      bool isSymmetric = true;
      fields.readRGrid(in, domain.unitCell(), isSymmetric);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }
};

TEST_BEGIN(WFieldContainerTest)
TEST_ADD(WFieldContainerTest, testAllocate_bcc)
TEST_ADD(WFieldContainerTest, testSetBasis_bcc)
TEST_ADD(WFieldContainerTest, testSetRGrid_1_bcc)
TEST_ADD(WFieldContainerTest, testSetRGrid_2_bcc)
TEST_ADD(WFieldContainerTest, testReadBasis_bcc)
TEST_ADD(WFieldContainerTest, testReadRGrid_1_bcc)
TEST_ADD(WFieldContainerTest, testReadRGrid_2_bcc)
TEST_END(WFieldContainerTest)

#endif
