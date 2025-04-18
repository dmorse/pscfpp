#ifndef POLYMER_STUB_TEST_H
#define POLYMER_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "PolymerStub.h"
#include <pscf/solvers/PolymerTmpl.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/VertexIterator.h>
#include <pscf/chem/EdgeIterator.h>
#include <util/containers/Pair.h>

#include <fstream>

using namespace Pscf;

class PolymerStubTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {  setVerbose(0); }
  
   void testReadParam(PolymerStub& p, std::string fileName) 
   {
      std::ifstream in;
      openInputFile(fileName, in);

      p.readParam(in);

      if (verbose() > 0) {
         std::cout << std::endl;
         p.writeParam(std::cout);
      }
    
      #if 0
      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << p.vertex(i).size() << "\n";
      }
      #endif

      if (verbose() > 0) {
         std::cout << "\nVertices: id, size, block ids\n";
         for (int i = 0; i < p.nVertex(); ++i) {
            std::cout << i << "  " << p.vertex(i).size();
            for (int j = 0; j < p.vertex(i).size(); ++j) {
               std::cout << "  " << p.vertex(i).inPropagatorId(j)[0];
            }
            std::cout << std::endl;
         }
      }

      #if 0
      for (int i = 0; i < p.nPropagator(); ++i) {
         std::cout << p.propagatorId(i)[0] << "  " 
                   << p.propagatorId(i)[1] << "\n";
      }
      #endif
     
      if (verbose() > 0) { 
         std::cout << "\nPropagator order:\n";
      }
      Pair<int> propId;
      PropagatorStub* propPtr = 0;
      for (int i = 0; i < p.nPropagator(); ++i) {
         propId = p.propagatorId(i);
         if (verbose() > 0) {
            std::cout << propId[0] << "  " << propId[1] << "\n";
         }
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

   int testVertexIterator(PolymerStub& p, int sourceId, int targetId) 
   {
      if (verbose() > 0) std::cout << "\nVertexIterator: \n";
      VertexIterator iter(p);
      iter.begin(sourceId, targetId);
      int nVertex = 1;
      while (iter.notEnd()) {
         if (verbose() > 0) {
            std::cout << iter.currentId() << std::endl;
         }
         ++iter;
         ++nVertex;
      }
      if (verbose() > 0) {
         std::cout << iter.currentId() << std::endl;
      }
      return nVertex;
   }

   int testEdgeIterator(PolymerStub& p, int sourceId, int targetId) 
   {
      if (verbose() > 0) std::cout << "\nEdgeIterator: \n";
      EdgeIterator iter(p);
      iter.begin(sourceId, targetId);
      int nEdge = 1;
      int edgeId, dirId, vertexId, vpId, vnId;
      while (iter.notEnd()) {
         edgeId = iter.currentEdgeId();
         dirId = iter.currentDirectionId();
         if (dirId == 0) {
            vpId = p.edge(edgeId).vertexId(0);
            vnId = p.edge(edgeId).vertexId(1);
         } else
         if (dirId == 1) {
            vpId = p.edge(edgeId).vertexId(1);
            vnId = p.edge(edgeId).vertexId(0);
         }
         if (nEdge > 1) {
            TEST_ASSERT(vpId == vertexId);
         }
         vertexId = iter.currentVertexId();
         TEST_ASSERT(vnId == vertexId);
         if (verbose() > 0) {
            std::cout << iter.currentEdgeId() << "  "
                      << iter.currentVertexId() << std::endl;
         }
         ++iter;
         ++nEdge;
      }
      if (verbose() > 0) {
         std::cout << iter.currentEdgeId() << "  " 
                   << iter.currentVertexId() << std::endl;
      }
      return nEdge;
   }

   // Unit test functions

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      PolymerStub p;
   } 

   void testReadParamDiblock() 
   {
      printMethod(TEST_FUNC);

      PolymerStub p;
      testReadParam(p, "in/PolymerDiblock");

      int nVertex = testVertexIterator(p, 2, 0);
      TEST_ASSERT(nVertex == 3);

      nVertex = testVertexIterator(p, 0, 2);
      TEST_ASSERT(nVertex == 3);

      int nEdge = testEdgeIterator(p, 0, 1); 
      TEST_ASSERT(nEdge == 2);
   }

   void testReadParamTriblock() 
   {
      printMethod(TEST_FUNC);

      PolymerStub p;
      testReadParam(p, "in/PolymerTriblock");

      int nVertex = testVertexIterator(p, 3, 0);
      TEST_ASSERT(nVertex == 4);

      int nEdge = testEdgeIterator(p, 2, 0); 
      TEST_ASSERT(nEdge == 3);
   }

   void testReadParamStar() 
   {
      printMethod(TEST_FUNC);

      /*
      * Star graph:
      *
      *           1
      *          /
      *     0 - 3
      *          \
      *           2
      *
      * Bond indices are the same as attached endpoints.
      */

      PolymerStub p;
      testReadParam(p, "in/PolymerStar");

      int nVertex = testVertexIterator(p, 2, 0);
      TEST_ASSERT(nVertex == 3);
 
      int nEdge = testEdgeIterator(p, 2, 1);
      TEST_ASSERT(nEdge == 2);
 
   }

   void testReadParamH() 
   {
      printMethod(TEST_FUNC);

      /*
      * H polymer graph topology (PolymerH file)
      *  
      *  0       2
      *   \     /
      *    4 - 5 
      *   /     \
      *  1       3
      * 
      * Bond indices are the same as indices of attached ends 
      * or the lesser of indices for the two attached vertices.
      */

      PolymerStub p;
      testReadParam(p, "in/PolymerH");

      int nVertex = testVertexIterator(p, 2, 1);
      TEST_ASSERT(nVertex == 4);
 
      int nEdge = testEdgeIterator(p, 2, 1); 
      TEST_ASSERT(nEdge == 3);
   }

};

TEST_BEGIN(PolymerStubTest)
TEST_ADD(PolymerStubTest, testConstructor)
TEST_ADD(PolymerStubTest, testReadParamDiblock)
TEST_ADD(PolymerStubTest, testReadParamTriblock)
TEST_ADD(PolymerStubTest, testReadParamStar)
TEST_ADD(PolymerStubTest, testReadParamH)
TEST_END(PolymerStubTest)

#endif
