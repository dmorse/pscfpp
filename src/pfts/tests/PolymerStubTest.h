#ifndef POLYMER_STUB_TEST_H
#define POLYMER_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "PolymerStub.h"
#include <util/containers/Pair.h>

#include <fstream>

using namespace Pfts;

class PolymerStubTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      PolymerStub p;
   } 

   void testReadParam() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/Polymer", in);

      PolymerStub p;
      p.readParam(in);

      p.writeParam(std::cout);

      std::cout << "\nVertices: id, size, block ids\n";
      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << i << "  " << p.vertex(i).size();
         for (int j = 0; j < p.vertex(i).size(); ++j) {
            std::cout << "  " << p.vertex(i).inPropagatorId(j)[0];
         }
         std::cout << "\n";
      }

      std::cout << "\nPropagator order:\n";
      Pair<int> propId;
      PropagatorStub* propPtr = 0;
      for (int i = 0; i < p.nPropagator(); ++i) {
         propId = p.propagatorId(i);
         std::cout << propId[0] << "  " << propId[1] << "\n";
         propPtr = &p.propagator(i);
         TEST_ASSERT(propPtr->block().id() == propId[0]);
         TEST_ASSERT(propPtr->directionId() == propId[1]);
         propPtr->setIsSolved(false);
      }

      // Check computation plan
      for (int i = 0; i < p.nPropagator(); ++i) {
         TEST_ASSERT(p.propagator(i).isReady());
         p.propagator(i).setIsSolved(true);
      }
 
   }

};

TEST_BEGIN(PolymerStubTest)
TEST_ADD(PolymerStubTest, testConstructor)
TEST_ADD(PolymerStubTest, testReadParam)
TEST_END(PolymerStubTest)

#endif
