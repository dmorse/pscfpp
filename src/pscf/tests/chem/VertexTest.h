#ifndef VERTEX_TEST_H
#define VERTEX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/BlockDescriptor.h>
#include <pscf/chem/Vertex.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class VertexTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Vertex v;
   } 

   void testAddBlock() {
      printMethod(TEST_FUNC);
      //printEndl();

      BlockDescriptor b;
      std::ifstream in;
      openInputFile("in/BlockDescriptor", in);

      in >> b;
      TEST_ASSERT(b.id() == 5);
      TEST_ASSERT(b.monomerId() == 0);
      TEST_ASSERT(b.vertexId(0) == 3);
      TEST_ASSERT(b.vertexId(1) == 4);
      TEST_ASSERT(eq(b.length(), 2.0));
      //std::cout << b << std::endl;

      Vertex v;
      v.setId(3);
      v.addBlock(b);
      TEST_ASSERT(v.size() == 1);
      TEST_ASSERT(v.outPropagatorId(0)[0] == 5);
      TEST_ASSERT(v.outPropagatorId(0)[1] == 0);
      TEST_ASSERT(v.inPropagatorId(0)[0] == 5);
      TEST_ASSERT(v.inPropagatorId(0)[1] == 1);
      //std::cout << v.inPropagatorId(0)[0] << "  "
      //          << v.inPropagatorId(0)[1] << "\n";
      //std::cout << v.outPropagatorId(0)[0] << "  "
      //          << v.outPropagatorId(0)[1] << "\n";
      
   }

};

TEST_BEGIN(VertexTest)
TEST_ADD(VertexTest, testConstructor)
TEST_ADD(VertexTest, testAddBlock)
TEST_END(VertexTest)

#endif
