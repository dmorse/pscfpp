#ifndef RPC_DOMAIN_TEST_H
#define RPC_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class DomainTest : public LogFileUnitTest 
{

   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
   }

   /*
   * Open and read file header to initialize Rp::Domain<D,CPT> system.
   */
   template <int D>
   void readHeader(std::string filename, Rp::Domain<D,CPT>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readRGridFieldHeader(in, nMonomer_);
      in.close();
   }

   #if 0
   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);

      std::ifstream in;
      openInputFile("in/Domain_old.prm", in);
      domain.readParam(in);
      in.close();

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.basis().nBasis() == 489);
   }
   #endif

   void testReadParamHeader() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);

      // Read parameter file
      std::ifstream in;
      openInputFile("in/Domain.prm", in);
      domain.readParam(in);
      in.close();

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.unitCell().lattice() == domain.lattice());
      TEST_ASSERT(domain.group().size() == 96);
      TEST_ASSERT(domain.basis().nBasis() == 0);

      // Read header and construct basis
      openInputFile("in/w_bcc.rf", in);
      bool isSymmetric;
      domain.fieldIo().readFieldHeader(in, nMonomer_, 
                                       domain.unitCell(), isSymmetric);
      in.close();

      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.group().size() > 1);
      TEST_ASSERT(domain.basis().nBasis() == 489);
   }

   void testReadHeader() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.basis().nBasis() == 489);
      TEST_ASSERT(nMonomer_ == 2);

      if (verbose() > 0) {
         // openLogFile("out/DomainTestReadHeader.log");
         Log::file() << "\n";
         Log::file() << "Cell  = " << domain.unitCell() << "\n";
         Log::file() << "Ngrid = " << domain.mesh().dimensions() << "\n";
         if (verbose() > 1) {
            domain.basis().outputStars(Log::file());
         }
      }

   }

   /*
   * Read the first cell parameter in the header of a file.
   */
   double readCellParam(std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      std::string word;
      double value = 0.0;
      while (in >> word) {
         if (word == "cell_param") {
            in >> value;
            break;
         }
      }
      in.close();
      return value;
   }

   /*
   * Write star and wave lists when the basis is constructed from the
   * header of a field file read into a local unit cell, as in the
   * field conversion commands, so that the unit cell of the Domain is
   * not set. The headers must give the cell of the field file.
   */
   void testWriteStarsWaves()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      fileMaster_.setRootPrefix(UnitTest::filePrefix());
      domain.setFileMaster(fileMaster_);
      std::ifstream in;
      openInputFile("in/Domain.prm", in);
      domain.readParam(in);
      in.close();

      {
         UnitCell<3> unitCell;
         bool isSymmetric;
         openInputFile("in/w_bcc.rf", in);
         domain.fieldIo().readFieldHeader(in, nMonomer_, unitCell,
                                          isSymmetric);
         in.close();
      }
      TEST_ASSERT(domain.basis().nBasis() == 489);
      TEST_ASSERT(!domain.unitCell().isInitialized());

      domain.writeStars("out/stars_bcc");
      domain.writeWaves("out/waves_bcc");
      double a = readCellParam("in/w_bcc.rf");
      TEST_ASSERT(a > 1.0);
      TEST_ASSERT(eq(readCellParam("out/stars_bcc"), a));
      TEST_ASSERT(eq(readCellParam("out/waves_bcc"), a));
   }

};

TEST_BEGIN(DomainTest)
//TEST_ADD(DomainTest, testReadParam)
TEST_ADD(DomainTest, testReadParamHeader)
TEST_ADD(DomainTest, testReadHeader)
TEST_ADD(DomainTest, testWriteStarsWaves)
TEST_END(DomainTest)

#endif
