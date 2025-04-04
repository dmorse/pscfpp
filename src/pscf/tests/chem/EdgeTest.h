#ifndef BLOCK_DESCRIPTOR_TEST_H
#define BLOCK_DESCRIPTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/Edge.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class EdgeTest : public UnitTest 
{

public:

   void setUp()
   {
      // setVerbose(1);
   }

   void tearDown()
   {
      PolymerModel::setModel(PolymerModel::Thread);
   }

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Edge v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);

      // Test Thread model
      PolymerModel::setModel(PolymerModel::Thread);

      std::ifstream in;
      openInputFile("in/Edge", in);

      Edge v;
      TEST_ASSERT(v.polymerType() == PolymerType::Branched);

      v.setId(5);
      in >> v;
      TEST_ASSERT(v.id() == 5);
      TEST_ASSERT(v.monomerId() == 0);
      TEST_ASSERT(v.vertexId(0) == 3);
      TEST_ASSERT(v.vertexId(1) == 4);
      TEST_ASSERT(eq(v.length(), 2.0));
      TEST_ASSERT(v.polymerType() == PolymerType::Branched);
      if (verbose() > 0) {
         printEndl();
         std::cout << v << std::endl ;
      }

      v.setPolymerType(PolymerType::Linear);
      v.setId(2);
      in >> v;
      TEST_ASSERT(v.id() == 2);
      TEST_ASSERT(v.monomerId() == 1);
      v.setVertexIds(2,3);
      TEST_ASSERT(v.vertexId(0) == 2);
      TEST_ASSERT(v.vertexId(1) == 3);
      TEST_ASSERT(eq(v.length(), 3.0));
      TEST_ASSERT(v.polymerType() == PolymerType::Linear);
      if (verbose() > 0) {
         std::cout << v << std::endl ;
      }

      v.setPolymerType(PolymerType::Branched);
      TEST_ASSERT(v.polymerType() == PolymerType::Branched);
      if (verbose() > 0) {
         std::cout << v << std::endl ;
      }

      // Test format for PolymerModel::Bead model
      PolymerModel::setModel(PolymerModel::Bead);

      // Linear polymer
      v.setPolymerType(PolymerType::Linear);
      v.setId(3);
      v.setVertexIds(3, 4);
      TEST_ASSERT(v.vertexId(0) == 3);
      TEST_ASSERT(v.vertexId(1) == 4);
      in >> v;
      TEST_ASSERT(v.id() == 3);
      TEST_ASSERT(v.monomerId() == 2);
      TEST_ASSERT(v.nBead() == 50);

      // Would cause run time error
      // cout << v.length();
      
   }

};

TEST_BEGIN(EdgeTest)
TEST_ADD(EdgeTest, testConstructor)
TEST_ADD(EdgeTest, testReadWrite)
TEST_END(EdgeTest)

#endif
