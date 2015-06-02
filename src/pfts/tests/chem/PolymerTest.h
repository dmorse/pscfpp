#ifndef POLYMER_TEST_H
#define POLYMER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pfts/chem/Polymer.h>
#include <pfts/chem/Vertex.h>
#include <pfts/chem/Block.h>

#include <fstream>

using namespace Pfts;
//using namespace Util;

class PolymerTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Polymer p;
   } 

   void testReadParam() {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/Polymer", in);

      Polymer p;
      p.readParam(in);
      p.writeParam(std::cout);

      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << p.vertex(i).size() << "\n";
      }
      
   }

};

TEST_BEGIN(PolymerTest)
TEST_ADD(PolymerTest, testConstructor)
TEST_ADD(PolymerTest, testReadParam)
TEST_END(PolymerTest)

#endif
