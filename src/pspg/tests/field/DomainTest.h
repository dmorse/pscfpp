#ifndef PSPG_DOMAIN_TEST_H
#define PSPG_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/field/Domain.h>
#include <pspg/field/FieldIo.h>
#include <pspg/field/FFT.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class DomainTest : public LogFileUnitTest 
{

   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   { setVerbose(0); }

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
      if (domain.basis().isInitialized()) {
         TEST_ASSERT(domain.basis().nBasis() == 489);
      }
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
         std::cout << "\n";
         std::cout << "Cell  = " << domain.unitCell() << "\n";
         std::cout << "Ngrid = " << domain.mesh().dimensions() << "\n";
         if (verbose() > 1) {
            domain.basis().outputStars(std::cout);
         }
      }

   }

};

TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testReadParam)
TEST_ADD(DomainTest, testReadHeader)
TEST_END(DomainTest)

#endif
