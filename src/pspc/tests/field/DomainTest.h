#ifndef PSPC_DOMAIN_TEST_H
#define PSPC_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <pscf/tests/LogFileUnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/Domain.h>
#include <pspc/field/FieldIo.h>
#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>
#include <pspc/field/FFT.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class DomainTest : public LogFileUnitTest 
{

   //std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
   }

   /*
   * Open and read file header to initialize Domain<D> system.
   */
   template <int D>
   void readHeader(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readFieldHeader(in, nMonomer_);
      in.close();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);

      std::ifstream in;
      openInputFile("in/Domain", in);
      domain.readParam(in);
      in.close();

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.basis().nBasis() == 489);
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

};

TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testReadParam)
TEST_ADD(DomainTest, testReadHeader)
TEST_END(DomainTest)

#endif
