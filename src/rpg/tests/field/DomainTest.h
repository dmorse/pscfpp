#ifndef RPG_DOMAIN_TEST_H
#define RPG_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>

#include <prdc/cuda/FFT.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Rpg;

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
      domain.readRGridFieldHeader(in, nMonomer_);
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
      TEST_ASSERT(domain.group().size() == 96);
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
      TEST_ASSERT(domain.group().size() == 96);
      TEST_ASSERT(domain.basis().isInitialized());
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
