#ifndef RPC_FIELD_IO_TEST_H
#define RPC_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>
#include <prdc/cpu/RFieldDftComparison.h>
#include <prdc/cpu/FFT.h>
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

class FieldIoTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      openLogFile("out/fieldIoTestLogFile");
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
   void allocateFields(int nMonomer, int nStar,
                            DArray< DArray<double> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(nStar);
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

   void testBasisIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_bcc.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_bcc.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_bcc.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }

   }

   void testBasisIo_c15_1() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_c15_1.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_c15_1.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_c15_1.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }

   }

   void testBasisIo_altG() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_altG.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_altG.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_altG.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_altG_fort() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_altG_fort.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_altG_fort.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_altG_fort.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      #if 1
      //setVerbose(1);
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 12) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(), 21, 12) << std::endl;
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-11);
      #endif

   }

   void testRGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      DArray< RField<3> > rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);

      readFields("in/w_bcc.rf", domain, rf_0);
      writeFields("out/w_bcc.rf", domain, rf_0);
      readFields("out/w_bcc.rf", domain, rf_1);

      RFieldComparison<3> comparison;
      comparison.compare(rf_0, rf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   }

   void testConvertBasisKGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }

   }

   void testConvertBasisRGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RField<3> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToRGrid(bf_0, rf_0);
      domain.fieldIo().convertRGridToBasis(rf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
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
      openOutputFile("out/stars_altG", out);
      domain.basis().outputStars(out);
      out.close();

      #if 1
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_altG.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      #endif

      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testConvertBasisKGridBasis_c15_1() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 2;

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      #if 0
      std::ofstream  out;
      openOutputFile("out/waves_c15_1", out);
      domain.basis().outputWaves(out);
      out.close();
      #endif

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_c15_1.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      // setVerbose(1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-10);
   }

   void testKGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      writeFields("out/w_bcc.kf", domain, kf_0);
      readFields("out/w_bcc.kf", domain, kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(kf_0, kf_1);
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

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      readFields("in/w_altG.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      writeFields("out/w_altG.kf", domain, kf_0);
      readFields("out/w_altG.kf", domain, kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(kf_0, kf_1);
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

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< RFieldDft<1> > kf_0, kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      readFields("in/w_lam.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      writeFields("out/w_lam.kf", domain, kf_0);
      readFields("out/w_lam.kf", domain, kf_1);

      RFieldDftComparison<1> comparison;
      comparison.compare(kf_0, kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
   }

   void testConvertBasisKGridRGridKGrid_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      DArray< RFieldDft<3> > kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_2);

      DArray< RField<3> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      kf_2 = kf_0;
      domain.fieldIo().convertKGridToRGrid(kf_0, rf_0);

      #if 0
      // Demonstrate that input kf_0 is modified by above (it is)
      RFieldDftComparison<3> check;
      check.compare(kf_2, kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(rf_0, kf_1);

      RFieldDftComparison<3> comparison;
      comparison.compare(kf_2, kf_1);
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

      // Read header
      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      // Allocate required fields
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      DArray< RFieldDft<3> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      DArray< RFieldDft<3> > kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_2);

      DArray< RField<3> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      // Read fields in basis format 
      readFields("in/w_c15_1.bf", domain, bf_0);

      // Convert basis -> kgrid -> rgrid -> kgrid -> basis
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      kf_2 = kf_0;
      domain.fieldIo().convertKGridToRGrid(kf_0, rf_0);

      #if 0
      // Demonstrate that input kf_0 is modified by above (it is)
      RFieldDftComparison<3> check;
      check.compare(kf_2, kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(rf_0, kf_1);

      RFieldDftComparison<3> Bcomparison;
      Bcomparison.compare(kf_2, kf_1);
      TEST_ASSERT(Bcomparison.maxDiff() < 1.0E-10);
      if (verbose() > 0) {
        std::cout  << std::endl;
        std::cout  << Dbl(Bcomparison.maxDiff(), 21, 13) << "\n";
        std::cout  << Dbl(Bcomparison.rmsDiff(), 21, 13) << "\n";
      }

      domain.fieldIo().convertKGridToBasis(kf_1, bf_1);
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

   }

   void testReplicate_bcc() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);
      IntVec<3> dimensions = domain.mesh().dimensions();

      DArray< RField<3> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("in/w_bcc.rf", domain, rf_0);

      IntVec<3> replicas;
      replicas[0] = 2;
      replicas[1] = 2;
      replicas[2] = 2;

      std::ofstream  out;
      openOutputFile("out/w_bcc_replica.rf", out);
      domain.fieldIo().replicateUnitCell(out, rf_0, 
                                         domain.unitCell(),
                                         replicas);

      // Read replicated field header 
      Domain<3> domain_rep;
      domain_rep.setFileMaster(fileMaster_);
      readHeader("out/w_bcc_replica.rf", domain_rep);
      IntVec<3> dimensions_rep = domain_rep.mesh().dimensions();
      for (int i=0; i < 3; ++i) {
         //TEST_ASSERT(domain_rep.mesh().dimension(i) 
         //            == 2*domain.mesh().dimension(i));
         TEST_ASSERT(dimensions_rep[i] == 2*dimensions[i]);
         TEST_ASSERT(dimensions[i] == 32);
      }

      // Read replicated field
      DArray< RField<3> > rf_1;
      allocateFields(nMonomer_, dimensions_rep, rf_1);
      readFields("out/w_bcc_replica.rf", domain_rep, rf_1);

      // Compare fields
      Mesh<3> mesh(dimensions);
      Mesh<3> mesh_rep(dimensions_rep);
      IntVec<3> p;
      IntVec<3> p_rep;
      int rank, rank_rep;
      for (p[0] = 0; p[0] < dimensions[0]; ++p[0]) {
         p_rep[0] = p[0] + dimensions[0];
         for (p[1] = 0; p[1] < dimensions[1]; ++p[1]) {
            p_rep[1] = p[1] + dimensions[1];
            for (p[2] = 0; p[2] < dimensions[2]; ++p[2]) {
               p_rep[2] = p[2] + dimensions[2];
               rank = mesh.rank(p);
               rank_rep = mesh_rep.rank(p_rep);
               //std::cout << Int(p[0]) << Int(p[1]) << Int(p[2])
               //          << Dbl(rf_0[0][rank])
               //          << Dbl(rf_1[0][rank_rep]) 
               //          << std::endl;
               for (int i=0; i < nMonomer_; ++i) {
                  TEST_ASSERT(eq( rf_0[i][rank], rf_1[i][rank_rep]));
               }
            }
         }
      }

   }

   void testExpand_lam_13() 
   {
      printMethod(TEST_FUNC);

      Domain<1> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);
      IntVec<1> dimensions = domain.mesh().dimensions();

      DArray< RField<1> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("in/w_lam.rf", domain, rf_0);

      int d = 3;
      DArray<int> newGridDimensions;
      newGridDimensions.allocate(2);
      newGridDimensions[0] = dimensions[0];
      newGridDimensions[1] = dimensions[0];

      std::ofstream  out;
      openOutputFile("out/w_lam_exp.rf", out);
      domain.fieldIo().expandRGridDimension(out, rf_0, 
                                  domain.unitCell(), d,
                                  newGridDimensions);

      // Read header for expanded field (d=3)
      Domain<3> domain_exp;
      domain_exp.setFileMaster(fileMaster_);
      readHeader("out/w_lam_exp.rf", domain_exp);
      IntVec<3> dimensions_exp = domain_exp.mesh().dimensions();
      dimensions_exp[0] = dimensions[0];
      dimensions_exp[1] = newGridDimensions[0];
      dimensions_exp[2] = newGridDimensions[1];

      // Read data for expanded field (d=3)
      DArray< RField<3> > rf_1;
      allocateFields(nMonomer_, dimensions_exp, rf_1);
      readFields("out/w_lam_exp.rf", domain_exp, rf_1);

      // Compare fields
      Mesh<3> mesh_exp(dimensions_exp);
      IntVec<3> p;
      int rank, rank_exp;
      for (p[0] = 0; p[0] < dimensions_exp[0]; ++p[0]) {
         rank = p[0];
         for (p[1] = 0; p[1] < dimensions_exp[1]; ++p[1]) {
            for (p[2] = 0; p[2] < dimensions_exp[2]; ++p[2]) {
               rank_exp = mesh_exp.rank(p);
               for (int i=0; i < nMonomer_; ++i) {
                  TEST_ASSERT(eq(rf_0[i][rank], rf_1[i][rank_exp]));
               }
            }
         }
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
TEST_ADD(FieldIoTest, testKGridIo_bcc)
TEST_ADD(FieldIoTest, testKGridIo_altG)
TEST_ADD(FieldIoTest, testKGridIo_lam)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_c15_1)
TEST_ADD(FieldIoTest, testReplicate_bcc)
TEST_ADD(FieldIoTest, testExpand_lam_13)
TEST_END(FieldIoTest)

#endif
