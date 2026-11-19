#ifndef RPC_W_FIELD_CONTAINER_TEST_H
#define RPC_W_FIELD_CONTAINER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>


#include <rp/field/WFields.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class WFieldsTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      openLogFile("out/wFieldContainerTest.log");
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

   // Open and read parameter header to initialize Rp::Domain<D,CPT> system.
   template <int D>
   void readParam(std::string filename, Rp::Domain<D,CPT>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readParam(in);
      in.close();
   }

   // Open and read file header to initialize Rp::Domain<D,CPT> system.
   template <int D>
   void readHeader(std::string filename, Rp::Domain<D,CPT>& domain)
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
                       DArray< RField<D,CPT> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(dimensions);
      }
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< DArray<double> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsBasis(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RField<D,CPT> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsRGrid(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< DArray<double> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsBasis(out, fields, domain.unitCell());
      out.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RField<D,CPT> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsRGrid(out, fields, domain.unitCell());
      out.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RFieldDft<D,CPT> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsKGrid(out, fields, domain.unitCell());
      out.close();
   }

   void testAllocate_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      Rp::WFields<3,CPT> fields;
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

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      //std::cout << "\n domain.basis().nBasis() = "
      //          << domain.basis().nBasis() << "\n";

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      Rp::WFields<3,CPT> fields;
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
      comparison.compare(bf, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testSetRGrid_1_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      readFields("in/w_bcc.rf", domain, rf);

      Rp::WFields<3,CPT> fields;
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

      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf, fields.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testSetRGrid_2_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      DArray< RField<3,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      domain.fieldIo().convertBasisToRGrid(bf, rf);

      Rp::WFields<3,CPT> fields;
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

   /*
   * Set fields in basis form in a container that holds fields in r-grid
   * form only, and has not allocated basis fields (READ_H_BASIS after
   * READ_H_RGRID, for the external fields).
   */
   void testBasisAfterRGrid_lam()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      fileMaster_.setRootPrefix(UnitTest::filePrefix());
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      readFields("in/w_lam.bf", domain, bf);
      DArray< RField<1,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      readFields("in/w_lam.rf", domain, rf);

      // readBasis after readRGrid
      Rp::WFields<1,CPT> fields;
      fields.setFieldIo(domain.fieldIo());
      fields.setNMonomer(nMonomer_);
      fields.setReadUnitCell(domain.unitCell());
      fields.readRGrid("in/w_lam.rf");
      TEST_ASSERT(!fields.isAllocatedBasis());
      bool thrown = false;
      try {
         fields.readBasis("in/w_lam.bf");
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(fields.isSymmetric());
      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      // setBasis after setRGrid with an asymmetric field
      Rp::WFields<1,CPT> fields2;
      fields2.setFieldIo(domain.fieldIo());
      fields2.setNMonomer(nMonomer_);
      fields2.setRGrid(rf);
      TEST_ASSERT(!fields2.isAllocatedBasis());
      thrown = false;
      try {
         fields2.setBasis(bf);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(fields2.isSymmetric());
      comparison.compare(bf, fields2.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testReadBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      Rp::WFields<3,CPT> fields;
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
      fields.setReadUnitCell(domain.unitCell());
      fields.readBasis(in);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      domain.fieldIo().convertRGridToBasis(fields.rgrid(), bf_1);
      comparison.compare(bf, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testReadRGrid_1_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      readFields("in/w_bcc.rf", domain, rf);

      Rp::WFields<3,CPT> fields;
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
      fields.setReadUnitCell(domain.unitCell());
      fields.readRGrid(in);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(!fields.isSymmetric());

      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf, fields.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testReadRGrid_2_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      DArray< RField<3,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      domain.fieldIo().convertBasisToRGrid(bf, rf);

      Rp::WFields<3,CPT> fields;
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
      fields.setReadUnitCell(domain.unitCell());
      fields.readRGrid(in, isSymmetric);
      TEST_ASSERT(fields.hasData());
      TEST_ASSERT(fields.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, fields.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }
};

TEST_BEGIN(WFieldsTest)
TEST_ADD(WFieldsTest, testAllocate_bcc)
TEST_ADD(WFieldsTest, testSetBasis_bcc)
TEST_ADD(WFieldsTest, testSetRGrid_1_bcc)
TEST_ADD(WFieldsTest, testSetRGrid_2_bcc)
TEST_ADD(WFieldsTest, testBasisAfterRGrid_lam)
TEST_ADD(WFieldsTest, testReadBasis_bcc)
TEST_ADD(WFieldsTest, testReadRGrid_1_bcc)
TEST_ADD(WFieldsTest, testReadRGrid_2_bcc)
TEST_END(WFieldsTest)

#endif
