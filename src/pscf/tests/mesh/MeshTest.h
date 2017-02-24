#ifndef PSCF_MESH_TEST_H
#define PSCF_MESH_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class MeshTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   void test3DMesh() {
      printMethod(TEST_FUNC);
      printEndl();

      Mesh<3> mesh;
      IntVec<3> d;
      d[0] = 2;
      d[1] = 1;
      d[2] = 3;
      mesh.setDimensions(d);

      IntVec<3> p;
      p[0] = 1;
      p[1] = 0;
      p[2] = 1;
      int rank = mesh.rank(p);
      IntVec<3> q = mesh.position(rank);
      TEST_ASSERT(q == p);

      //std::cout << "position = " << p << std::endl;
      //std::cout << "rank     = " << rank << std::endl;
      //std::cout << "position = " << q << std::endl;

   }

};

TEST_BEGIN(MeshTest)
TEST_ADD(MeshTest, test3DMesh)
TEST_END(MeshTest)

#endif
