#ifndef BLOCK_DESCRIPTOR_TEST_H
#define BLOCK_DESCRIPTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/Edge.h>

#include <cmath>
#include <fstream>
#include <sstream>

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

   /*
   * A length written by operator << and read by operator >> must be
   * the same double. Lengths were written with 13 significant digits.
   */
   void testWriteReadLength()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      double lengths[3];
      lengths[0] = 1.0/3.0;
      lengths[1] = 0.2 + 0.085;   // 0.28500000000000003
      lengths[2] = 2.0*std::acos(-1.0);

      Edge v, u;
      v.setMonomerId(1);
      v.setVertexIds(2, 5);
      for (int i = 0; i < 3; ++i) {
         v.setLength(lengths[i]);
         std::stringstream stream;
         stream << v;
         stream >> u;
         TEST_ASSERT(u.monomerId() == 1);
         TEST_ASSERT(u.vertexId(0) == 2);
         TEST_ASSERT(u.vertexId(1) == 5);
         TEST_ASSERT(u.length() == lengths[i]);
      }

      // 15 significant digits if these give the same double
      v.setLength(0.8);
      std::stringstream stream;
      stream << v;
      TEST_ASSERT(stream.str().find("8.00000000000000e-01")
                  != std::string::npos);
   }

   /*
   * Output to a stream set to fixed format: setf(scientific) without
   * the floatfield mask gave hexadecimal output, which operator >>
   * cannot read.
   */
   void testWriteFixedStream()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      Edge v, u;
      v.setMonomerId(1);
      v.setVertexIds(2, 5);
      v.setLength(1.25);
      std::stringstream stream;
      stream << std::fixed << v;
      TEST_ASSERT(stream.str().find("0x") == std::string::npos);
      stream >> u;
      TEST_ASSERT(u.length() == 1.25);
      TEST_ASSERT(u.vertexId(1) == 5);
   }

};

TEST_BEGIN(EdgeTest)
TEST_ADD(EdgeTest, testConstructor)
TEST_ADD(EdgeTest, testReadWrite)
TEST_ADD(EdgeTest, testWriteReadLength)
TEST_ADD(EdgeTest, testWriteFixedStream)
TEST_END(EdgeTest)

#endif
