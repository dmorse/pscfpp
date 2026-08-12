#ifndef CPC_DOMAIN_TEST_H
#define CPC_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cpc/field/Domain.h>
#include <cpc/field/FieldIo.h>

#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldDft.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Cpc;

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
   }

   void testReadParamHeader() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
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

      // Read field file header
      openInputFile("in/w_bcc.rf", in);
      domain.fieldIo().readFieldHeader(in, nMonomer_, domain.unitCell());
      in.close();

      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.lattice() == UnitCell<3>::Cubic);
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
      TEST_ASSERT(nMonomer_ == 2);

      //setVerbose(1);
      if (verbose() > 0) {
         openLogFile("out/DomainTestReadHeader.log");
         Log::file() << "\n";
         Log::file() << "Cell  = " << domain.unitCell() << "\n";
         Log::file() << "Ngrid = " << domain.mesh().dimensions() << "\n";
      }

   }

};

TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testReadParam)
TEST_ADD(DomainTest, testReadParamHeader)
TEST_ADD(DomainTest, testReadHeader)
TEST_END(DomainTest)

#endif
