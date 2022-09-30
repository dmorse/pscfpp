#ifndef PSPC_MASK_TEST_H
#define PSPC_MASK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/RFieldComparison.h>

#include <pspc/field/Mask.h>
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

class MaskTest : public UnitTest 
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
      int nBasis = domain.basis().nBasis();
      IntVec<3> dimensions = domain.mesh().dimensions();

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, nBasis, bf);
      readFields("in/w_bcc.bf", domain, bf);

      Mask<3> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocate(nBasis, dimensions);
      TEST_ASSERT(mask.isAllocated());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      mask.setBasis(bf[0]);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf[0], mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray<double> bf_1;
      bf_1.allocate(nBasis);
      domain.fieldIo().convertRGridToBasis(mask.rgrid(), bf_1);
      comparison.compare(bf_1, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testSetRGrid_1_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<3> dimensions = domain.mesh().dimensions();

      DArray< RField<3> > rf;
      allocateFields(nMonomer_, dimensions, rf);
      TEST_ASSERT(rf.capacity() == nMonomer_);
      readFields("in/w_bcc.rf", domain, rf);

      Mask<3> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocate(nBasis, dimensions);
      TEST_ASSERT(mask.isAllocated());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      mask.setRGrid(rf[0]);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      RFieldComparison<3> comparison;
      comparison.compare(rf[0], mask.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testSetRGrid_2_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<3> dimensions = domain.mesh().dimensions();

      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, nBasis, bf);
      TEST_ASSERT(bf.capacity() == nMonomer_);
      readFields("in/w_bcc.bf", domain, bf);

      RField<3> rf;
      rf.allocate(dimensions);
      domain.fieldIo().convertBasisToRGrid(bf[0], rf);

      Mask<3> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocate(nBasis, dimensions);
      TEST_ASSERT(mask.isAllocated());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      bool isSymmetric = true;
      mask.setRGrid(rf, isSymmetric);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf[0], mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

};

TEST_BEGIN(MaskTest)
TEST_ADD(MaskTest, testSetBasis_bcc)
TEST_ADD(MaskTest, testSetRGrid_1_bcc)
TEST_ADD(MaskTest, testSetRGrid_2_bcc)
TEST_END(MaskTest)

#endif
