#ifndef RPG_FIELD_IO_TEST_H
#define RPG_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldDft.h>
#include <prdc/cuda/RFieldComparison.h>
#include <prdc/cuda/RFieldDftComparison.h>
#include <prdc/cuda/FFT.h>

#include <prdc/crystal/BFieldComparison.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/inter/Interaction.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class FieldIoTest : public LogFileUnitTest 
{

   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      fileMaster_.setRootPrefix(UnitTest::filePrefix());
   }

   /*
   * Open and read parameter header to initialize Domain<D> system.
   */
   template <int D>
   void readParam(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readParam(in);
      in.close();
   }

   /*
   * Open and read file header to initialize Domain<D> system.
   */
   template <int D>
   void readHeader(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readRGridFieldHeader(in, nMonomer_);
      in.close();
   }

   // Allocate an array of fields in symmetry adapated format
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

   // Allocate an array of k-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
                            DArray< RFieldDft<D> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(dimensions);
      }
   }

   template <int D>
   void readFieldsBasis(std::string filename, Domain<D>& domain,
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
   void readFields(std::string filename, Domain<D>& domain,
                   DArray< RFieldDft<D> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsKGrid(in, fields, domain.unitCell());
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

   void testReadHeader() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.basis().nBasis() == 489);
      //TEST_ASSERT(nMonomer_ == 2);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell  = " << domain.unitCell() << "\n";
         std::cout << "Ngrid = " << domain.mesh().dimensions() << "\n";
         if (verbose() > 1) {
            domain.basis().outputStars(std::cout);
         }
      }

      DArray< DArray<double> > fb;
      allocateFields(nMonomer_, domain.basis().nBasis(), fb);

      DArray< RField<3> >  fr;
      allocateFields(nMonomer_, domain.mesh().dimensions(), fr);

      DArray< RFieldDft<3> > fk;
      allocateFields(nMonomer_, domain.mesh().dimensions(), fk);

   }

   void testBasisIo3D(std::string rf, std::string bf)
   {
      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/" + rf, domain);

      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, nBasis, d_bf_0);      

      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, nBasis, d_bf_1);

      std::ifstream in;
      openInputFile("in/" + bf, in);
      domain.fieldIo().readFieldsBasis(in, d_bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/" + bf, out);
      domain.fieldIo().writeFieldsBasis(out, d_bf_0, domain.unitCell());
      out.close();

      openInputFile("out/" + bf, in);
      domain.fieldIo().readFieldsBasis(in, d_bf_1, domain.unitCell());
      in.close();

      // Perform comparison
      BFieldComparison comparison;
      comparison.compare(d_bf_0, d_bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }

   }

   void testBasisIo_bcc() 
   {
      printMethod(TEST_FUNC);

      testBasisIo3D("w_bcc.rf", "w_bcc.bf");
   }

   void testBasisIo_c15_1() 
   {
      printMethod(TEST_FUNC);

      testBasisIo3D("c_c15_1.rf","w_c15_1.bf");
   }

   void testBasisIo_altG() 
   {
      printMethod(TEST_FUNC);

      testBasisIo3D("w_altG.rf", "w_altG.bf");
   }

   void testBasisIo_altG_fort() 
   {
      printMethod(TEST_FUNC);

      testBasisIo3D("w_altG.rf", "w_altG_fort.bf");
   }

   void testRGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3> > d_rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_rf_0);
      DArray< RField<3> > d_rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_rf_1);

      readFields("in/w_bcc.rf", domain, d_rf_0);
      writeFields("out/w_bcc.rf", domain, d_rf_0);
      readFields("out/w_bcc.rf", domain, d_rf_1);

      RFieldComparison<3> comparison;
      comparison.compare(d_rf_0, d_rf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
   }

   void testConvertBasisKGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, nBasis, d_bf_0);
      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, nBasis, d_bf_1);
      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      readFieldsBasis("in/w_bcc.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);
      domain.fieldIo().convertKGridToBasis(d_kf_0, d_bf_1);

      std::ofstream  out;
      openOutputFile("out/w_bcc_convert.bf", out);
      domain.fieldIo().writeFieldsBasis(out, d_bf_1, domain.unitCell());
      out.close();

      BFieldComparison comparison;
      comparison.compare(d_bf_0, d_bf_1);
      //setVerbose(1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testConvertBasisRGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, nBasis, d_bf_0);
      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, nBasis, d_bf_1);
      DArray< RField<3> > d_rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_rf_0);

      readFieldsBasis("in/w_bcc.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToRGrid(d_bf_0, d_rf_0);
      domain.fieldIo().convertRGridToBasis(d_rf_0, d_bf_1);

      BFieldComparison comparison;
      comparison.compare(d_bf_0, d_bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
   }

   void testConvertBasisKGridBasis_altG() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 3;

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      std::ofstream  out;
      openOutputFile("out/altG.stars", out);
      domain.basis().outputStars(out);
      out.close();

      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, nBasis, d_bf_0);
      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, nBasis, d_bf_1);
      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      readFieldsBasis("in/w_altG.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);
      domain.fieldIo().convertKGridToBasis(d_kf_0, d_bf_1);

      openOutputFile("out/w_altG_convert.bf", out);
      domain.fieldIo().writeFieldsBasis(out, d_bf_1, domain.unitCell());
      out.close();

      BFieldComparison comparison;
      comparison.compare(d_bf_0, d_bf_1);

      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }

   }

   void testConvertBasisKGridBasis_c15_1() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, nBasis, d_bf_0);
      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, nBasis, d_bf_1);
      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      readFieldsBasis("in/w_c15_1.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);
      domain.fieldIo().convertKGridToBasis(d_kf_0, d_bf_1);

      std::ofstream  out;
      openOutputFile("out/w_c15_1_convert.bf", out);
      domain.fieldIo().writeFieldsBasis(out, d_bf_1, domain.unitCell());
      out.close();

      BFieldComparison comparison;
      comparison.compare(d_bf_0, d_bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }

   }

   void testEstimateWBasis() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 2;

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      domain.fieldIo().setNMonomer(nMonomer_);
      readHeader("in/c_c15_1.rf", domain);

      Interaction v;
      v.setNMonomer(nMonomer_);
      std::ifstream in;
      openInputFile("in/Interaction", in);
      v.readParam(in);
      in.close();

      domain.fieldIo().estimateWBasis("in/c_c15_1.bf","out/w_c15_est.bf",
                                      v.chi());

      DArray< DArray<double> > bf_0, bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      readFieldsBasis("out/w_c15_est.bf", domain, bf_0);
      readFieldsBasis("in/w_c15_est.bf", domain, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-8);
   }

   void testKGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_0);

      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      DArray< RFieldDft<3> > d_kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_1);

      readFieldsBasis("in/w_bcc.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);

      writeFields("out/w_bcc.kf", domain, d_kf_0);
      readFields("out/w_bcc.kf", domain, d_kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(d_kf_0, d_kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-11);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
   }

   void testKGridIo_altG() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_0);

      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      DArray< RFieldDft<3> > d_kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_1);

      readFieldsBasis("in/w_altG.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);

      writeFields("out/w_altG.kf", domain, d_kf_0);
      readFields("out/w_altG.kf", domain, d_kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(d_kf_0, d_kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-11);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
   }

   void testKGridIo_lam() 
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_0);

      DArray< RFieldDft<1> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      DArray< RFieldDft<1> > d_kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_1);

      readFieldsBasis("in/w_lam.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);

      writeFields("out/w_lam.kf", domain, d_kf_0);
      readFields("out/w_lam.kf", domain, d_kf_1);

      RFieldDftComparison<1> comparison;
      comparison.compare(d_kf_0, d_kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
   }

   void testConvertBasisKGridRGridKGrid_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_0);

      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_1);

      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      DArray< RFieldDft<3> > d_kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_1);

      DArray< RFieldDft<3> > d_kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_2);

      DArray< RField<3> > d_rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_rf_0);

      readFieldsBasis("in/w_bcc.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);

      d_kf_2 = d_kf_0;
      domain.fieldIo().convertKGridToRGrid(d_kf_0, d_rf_0);

      #if 0
      // Demonstrate that input d_kf_0 is destroyed/overwritten by above
      RFieldDftComparison<3> check;
      check.compare(d_kf_2, d_kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(d_rf_0, d_kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(d_kf_2, d_kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
      if (verbose() > 0) {
        std::cout  << std::endl;
        std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
        std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
   }

   void testConvertBasisKGridRGridKGrid_c15_1() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      DArray< DArray<double> > d_bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_0);

      DArray< DArray<double> > d_bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), d_bf_1);

      DArray< RFieldDft<3> > d_kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_0);

      DArray< RFieldDft<3> > d_kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_1);

      DArray< RFieldDft<3> > d_kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_kf_2);

      DArray< RField<3> > d_rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), d_rf_0);

      readFieldsBasis("in/w_c15_1.bf", domain, d_bf_0);
      domain.fieldIo().convertBasisToKGrid(d_bf_0, d_kf_0);

      d_kf_2 = d_kf_0;
      domain.fieldIo().convertKGridToRGrid(d_kf_0, d_rf_0);

      #if 0
      // Demonstrate that input d_kf_0 is destroyed/overwritten by above
      RFieldDftComparison<3> check;
      check.compare(d_kf_2, d_kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(d_rf_0, d_kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(d_kf_2, d_kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
      if (verbose() > 0) {
        std::cout  << std::endl;
        std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
        std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
   }

};

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testReadHeader)
TEST_ADD(FieldIoTest, testBasisIo_bcc)
TEST_ADD(FieldIoTest, testBasisIo_c15_1)
TEST_ADD(FieldIoTest, testBasisIo_altG)
TEST_ADD(FieldIoTest, testBasisIo_altG_fort)
TEST_ADD(FieldIoTest, testRGridIo_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_bcc)
TEST_ADD(FieldIoTest, testConvertBasisRGridBasis_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_altG) 
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_c15_1)
TEST_ADD(FieldIoTest, testEstimateWBasis)
TEST_ADD(FieldIoTest, testKGridIo_bcc)
TEST_ADD(FieldIoTest, testKGridIo_altG)
TEST_ADD(FieldIoTest, testKGridIo_lam)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_c15_1)
TEST_END(FieldIoTest)

#endif
