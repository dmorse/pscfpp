#ifndef RPG_MASK_TEST_H
#define RPG_MASK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/field/Mask.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>

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
using namespace Pscf::Rpg;

class MaskTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;

public:

   void setUp()
   {
      setVerbose(0);
      openLogFile("out/maskTestLogFile");
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
      int nm;
      domain.readRGridFieldHeader(in, nm);
      in.close();
   }

   template <int D>
   void readField(std::string filename, Domain<D>& domain,
                   DArray<double>& field)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldBasis(in, field, domain.unitCell());
      in.close();
   }

   template <int D>
   void readField(std::string filename, Domain<D>& domain,
                  RField<D>& field)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldRGrid(in, field, domain.unitCell());
      in.close();
   }

   void testSetBasis() 
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      DArray<double> bf;
      bf.allocate(nBasis);
      readField("in/mask.bf", domain, bf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      mask.setBasis(bf);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray<double> bf_1;
      bf_1.allocate(nBasis);
      domain.fieldIo().convertRGridToBasis(mask.rgrid(), bf_1);
      comparison.compare(bf_1, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testSetRGrid_1() 
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      RField<1> rf;
      rf.allocate(dimensions);
      TEST_ASSERT(rf.capacity() == domain.mesh().size());
      readField("in/mask.rf", domain, rf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      mask.setRGrid(rf);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      RFieldComparison<1> comparison;
      comparison.compare(rf, mask.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testSetRGrid_2() 
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      DArray<double> bf;
      bf.allocate(nBasis);
      readField("in/mask.bf", domain, bf);

      RField<1> rf;
      rf.allocate(dimensions);
      domain.fieldIo().convertBasisToRGrid(bf, rf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      bool isSymmetric = true;
      mask.setRGrid(rf, isSymmetric);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testReadBasis()
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      DArray<double> bf;
      bf.allocate(nBasis);
      readField("in/mask.bf", domain, bf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      std::ifstream in;
      openInputFile("in/mask.bf", in);
      mask.readBasis(in, domain.unitCell());
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      BFieldComparison comparison;
      comparison.compare(bf, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      DArray<double> bf_1;
      bf_1.allocate(nBasis);
      domain.fieldIo().convertRGridToBasis(mask.rgrid(), bf_1);
      comparison.compare(bf_1, mask.basis());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testReadRGrid_1()
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      RField<1> rf;
      rf.allocate(dimensions);
      TEST_ASSERT(rf.capacity() == domain.mesh().size());
      readField("in/mask.rf", domain, rf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      std::ifstream in;
      openInputFile("in/mask.rf", in);
      mask.readRGrid(in, domain.unitCell());
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      RFieldComparison<1> comparison;
      comparison.compare(rf, mask.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testReadRGrid_2()
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      RField<1> rf;
      rf.allocate(dimensions);
      TEST_ASSERT(rf.capacity() == domain.mesh().size());
      readField("in/mask.rf", domain, rf);

      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());

      std::ifstream in;
      openInputFile("in/mask.rf", in);
      mask.readRGrid(in, domain.unitCell(), true);
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());

      RFieldComparison<1> comparison;
      comparison.compare(rf, mask.rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
   }

   void testPhiTot()
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/mask.rf", domain);
      int nBasis = domain.basis().nBasis();
      IntVec<1> dimensions = domain.mesh().dimensions();

      // Create empty mask object, check phiTot
      Mask<1> mask;
      mask.setFieldIo(domain.fieldIo());
      mask.allocateBasis(nBasis);
      mask.allocateRGrid(dimensions);
      TEST_ASSERT(mask.isAllocatedBasis());
      TEST_ASSERT(mask.isAllocatedRGrid());
      TEST_ASSERT(!mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());
      TEST_ASSERT(eq(mask.phiTot(), 1.0));

      // Read unsymmetrized r-grid, check phiTot
      std::ifstream in;
      openInputFile("in/mask.rf", in);
      mask.readRGrid(in, domain.unitCell());
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(!mask.isSymmetric());
      TEST_ASSERT(eq(mask.phiTot(), 8.9461021637e-01));

      // Read basis, check phiTot
      std::ifstream in2;
      openInputFile("in/mask.bf", in2);
      mask.readBasis(in2, domain.unitCell());
      TEST_ASSERT(mask.hasData());
      TEST_ASSERT(mask.isSymmetric());
      TEST_ASSERT(eq(mask.phiTot(), mask.basis()[0]));
      TEST_ASSERT(eq(mask.phiTot(), 8.9461021637e-01));
   }

};

TEST_BEGIN(MaskTest)
TEST_ADD(MaskTest, testSetBasis)
TEST_ADD(MaskTest, testSetRGrid_1)
TEST_ADD(MaskTest, testSetRGrid_2)
TEST_ADD(MaskTest, testReadBasis)
TEST_ADD(MaskTest, testReadRGrid_1)
TEST_ADD(MaskTest, testReadRGrid_2)
TEST_ADD(MaskTest, testPhiTot)
TEST_END(MaskTest)

#endif
