#ifndef POLYMER_DESCRIPTOR_TEST_H
#define POLYMER_DESCRIPTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/PolymerDescriptor.h>
#include <pscf/Vertex.h>
#include <pscf/Block.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class PolymerDescriptorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      PolymerDescriptor p;
   } 

   void testReadParam() {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/PolymerDescriptor", in);

      PolymerDescriptor p;
      p.readParam(in);
      p.writeParam(std::cout);

      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << p.vertex(i).size() << "\n";
      }

      for (int i = 0; i < p.nPropagator(); ++i) {
         std::cout << p.propagatorId(i)[0] << "  " 
                   << p.propagatorId(i)[1] << "\n";
      }
      
   }

};

TEST_BEGIN(PolymerDescriptorTest)
TEST_ADD(PolymerDescriptorTest, testConstructor)
TEST_ADD(PolymerDescriptorTest, testReadParam)
TEST_END(PolymerDescriptorTest)

#endif
