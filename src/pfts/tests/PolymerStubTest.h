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

      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << p.vertex(i).size() << "\n";
      }

      Pair<int> propId;
      PropagatorStub* propPtr = 0;
      for (int i = 0; i < p.nPropagator(); ++i) {
         propId = p.propagatorId(i);
         propPtr = &p.propagator(propId[0], propId[1]);
         propPtr->setIsComplete(false);
      }

      for (int i = 0; i < p.nPropagator(); ++i) {
         propId = p.propagatorId(i);
         std::cout << propId[0] << "  " << propId[1] << "\n";
         propPtr = &p.propagator(propId[0], propId[1]);
         propPtr->setIsComplete(false);
      }
      
   }

};

TEST_BEGIN(PolymerStubTest)
TEST_ADD(PolymerStubTest, testConstructor)
TEST_ADD(PolymerStubTest, testReadParam)
TEST_END(PolymerStubTest)

#endif
