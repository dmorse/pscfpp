#ifndef PSCF_MESH_ITERATOR_FORTRAN_TEST_H
#define PSCF_MESH_ITERATOR_FORTRAN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class MeshIteratorFortranTest : public UnitTest 
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
      MeshIteratorFortran<3> iter;
      iter.setDimensions(d);

      IntVec<3> x;  // current position
      IntVec<3> xp; // previous position
      int i ;       // current rank
      int r;
      
      if (verbose() > 0) {
         std::cout << iter.size() << std::endl;
         std::cout << d << std::endl;
         std::cout << iter.offsets() << std::endl;
         std::cout << std::endl;
      }
      int j = 0; // Sum of rank values
      for (iter.begin(); !iter.atEnd(); ++iter) {
         if (iter.rank() > 0) {
            xp = x;
         }
         x = iter.position();
         i = iter.rank();
         j += i;
         r = x[0]*d[1]*d[2] + x[1]*d[2] + x[2];
         if (verbose() > 0) {
            std::cout << x << "            " << i << std::endl;
         }
         TEST_ASSERT(i == r);
         if (i == 0) {
            for (int k = 0; k < 3; ++k) {
               TEST_ASSERT(x[k] == 0);
            }
         } else {
            if (x[0] != 0) {
               TEST_ASSERT(xp[0] == x[0] - 1);
               TEST_ASSERT(xp[1] == x[1]);
               TEST_ASSERT(xp[2] == x[2]);
            } else {
               if (x[1] != 0) {
                  TEST_ASSERT(xp[1] == x[1] - 1);
                  TEST_ASSERT(xp[2] == x[2]);
               } else {
                  TEST_ASSERT(xp[2] == x[2] - 1);
               }
            }
         }
      }

      // Each rank should appear once in range 0, ..., s-1
      // The sum of integers from 0, ..., s-1 should be s(s-1)/2
      int s = iter.size();
      TEST_ASSERT(j == s*(s-1)/2);

   }

   void test2D(IntVec<2>& d) 
   {
      MeshIteratorFortran<2> iter;
      iter.setDimensions(d);

      IntVec<2> x;  // current position
      IntVec<2> xp; // previous position
      int i ;       // current rank
      int r;        // manually computed rank
      if (verbose() > 0) {
         std::cout << iter.size() << std::endl;
         std::cout << d << std::endl;
         std::cout << iter.offsets() << std::endl;
         std::cout << std::endl;
      }
      int j = 0; // Sum of rank values
      for (iter.begin(); !iter.atEnd(); ++iter) {
         if (iter.rank() > 0) {
            xp = x;
         }
         x = iter.position();
         i = iter.rank();
         j += i;
         r = x[0]*d[1] + x[1];
         if (verbose() > 0) {
            std::cout << x << "            " << i << std::endl;
         }
         TEST_ASSERT(r == i);
         if (iter.rank() == 0) {
            TEST_ASSERT(x[1] == 0);
            TEST_ASSERT(x[0] == 0);
            TEST_ASSERT(i == 0);
         } else {
            if (x[0] != 0) {
               TEST_ASSERT(xp[0] == x[0] - 1);
               TEST_ASSERT(xp[1] == x[1]);
            } else {
               TEST_ASSERT(xp[1] == x[1] - 1);
            }
         }
      }

      // Each rank should appear once in range 0, ..., s-1
      // The sum of integers from 0, ..., s-1 should be s(s-1)/2
      int s = iter.size();
      TEST_ASSERT(j == s*(s-1)/2);
   }

   void testMeshIteratorFortran3D() 
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

   void testMeshIteratorFortran2D() 
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

TEST_BEGIN(MeshIteratorFortranTest)
TEST_ADD(MeshIteratorFortranTest, testMeshIteratorFortran3D)
TEST_ADD(MeshIteratorFortranTest, testMeshIteratorFortran2D)
TEST_END(MeshIteratorFortranTest)

#endif
