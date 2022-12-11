#ifndef PSCF_MESH_ITERATOR_TEST_H
#define PSCF_MESH_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class MeshIteratorTest : public UnitTest 
{

public:

   void setUp()
   {
      //setVerbose(1);
   }

   void tearDown()
   {}
 
   void test3D(IntVec<3>& d) 
   {
      MeshIterator<3> iter;
      iter.setDimensions(d);

      IntVec<3> x;  // current position
      IntVec<3> xp; // previous position
      int i ;       // current rank
      int ip ;      // previous rank
      for (iter.begin(); !iter.atEnd(); ++iter) {
         if (verbose() > 0) {
            std::cout << iter.rank() << "   " 
                        << iter.position() << std::endl;
         }
         if (iter.rank() == 0) {
            x = iter.position();
            i = iter.rank();
            for (int k = 0; k < 3; ++k) {
               TEST_ASSERT(x[k] == 0);
            }
         } else {
            xp = x;
            ip = i;
            x = iter.position();
            i = iter.rank();
            TEST_ASSERT(ip == i - 1);
            if (x[2] != 0) {
               TEST_ASSERT(xp[2] == x[2] - 1);
               TEST_ASSERT(xp[1] == x[1]);
               TEST_ASSERT(xp[0] == x[0]);
            } else {
               TEST_ASSERT(xp[2] == d[2] - 1);
               if (x[1] != 0) {
                  TEST_ASSERT(xp[1] == x[1] - 1);
               } else {
                  TEST_ASSERT(xp[1] == d[1] - 1);
                  TEST_ASSERT(xp[0] == x[0] - 1);
               }
            }
         }
      }

   }

   void test2D(IntVec<2>& d) 
   {
      MeshIterator<2> iter;
      iter.setDimensions(d);

      IntVec<2> x;  // current position
      IntVec<2> xp; // previous position
      int i ;       // current rank
      int ip ;      // previous rank
      for (iter.begin(); !iter.atEnd(); ++iter) {
         if (verbose() > 0) {
            std::cout << iter.rank() << "   " 
                      << iter.position() << std::endl;
         }
         if (iter.rank() == 0) {
            x = iter.position();
            i = iter.rank();
            TEST_ASSERT(x[1] == 0);
            TEST_ASSERT(x[0] == 0);
         } else {
            xp = x;
            ip = i;
            x = iter.position();
            i = iter.rank();
            TEST_ASSERT(ip == i - 1);
            if (x[1] != 0) {
               TEST_ASSERT(xp[1] == x[1] - 1);
               TEST_ASSERT(xp[0] == x[0]);
            } else {
               TEST_ASSERT(xp[1] == d[1] - 1);
               TEST_ASSERT(xp[0] == x[0] - 1);
            }
         }
      }

   }

   void testMeshIterator3D() 
   {
      printMethod(TEST_FUNC);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      IntVec<3> d;
      d[0] = 2;
      d[1] = 1;
      d[2] = 3;
      test3D(d);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      d[0] = 1;
      d[1] = 3;
      d[2] = 2;
      test3D(d);
   }

   void testMeshIterator2D() 
   {
      printMethod(TEST_FUNC);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      IntVec<2> twoD;
      twoD[0] = 2;
      twoD[1] = 3;
      test2D(twoD);
   }

};

TEST_BEGIN(MeshIteratorTest)
TEST_ADD(MeshIteratorTest, testMeshIterator3D)
TEST_ADD(MeshIteratorTest, testMeshIterator2D)
TEST_END(MeshIteratorTest)

#endif
