#ifndef PSSP_BASIS_TEST_H
#define PSSP_BASIS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

#include <pssp/basis/Basis.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class BasisTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor()
   {
      Basis<3> sampleBasis;
      TEST_ASSERT(eq(sampleBasis.nWave(),0));
      TEST_ASSERT(eq(sampleBasis.nStar(),0));
      TEST_ASSERT(eq(sampleBasis.nBasis(),0));
   }

   void testMakeBasis(){
      Basis<3> sampleBasis;

      UnitCell<3> unitCell;
      //make unitcell
      std::ifstream in;
      openInputFile("in/UnitCell", in);
      in >> unitCell;
      IntVec<3> d;
      in >> d;
      in.close();

      //make mesh object
      Mesh<3> mesh(d);
      std::string spaceGroup = "I";
      sampleBasis.makeBasis(mesh, unitCell, spaceGroup);

      //watch as the world burns
   }

};

TEST_BEGIN(BasisTest)
TEST_ADD(BasisTest, testConstructor)
TEST_ADD(BasisTest, testMakeBasis)
TEST_END(BasisTest)

#endif
