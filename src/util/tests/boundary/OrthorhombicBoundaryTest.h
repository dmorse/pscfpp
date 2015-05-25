#ifndef ORTHORHOMBIC_BOUNDARY_TEST
#define ORTHORHOMBIC_BOUNDARY_TEST

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/boundary/OrthorhombicBoundary.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/random/Random.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <fstream>

using namespace Util;

class OrthorhombicBoundaryTest : public UnitTest 
{

private:

   OrthorhombicBoundary boundary;

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testInitialize() 
   {
      printMethod(TEST_FUNC);

      Vector L, Lp;

      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      boundary.setOrthorhombic(L);
      Lp = boundary.lengths();


      // Assertions
      TEST_ASSERT(boundary.isValid());

      #if 0
      // Verbose output
      if (verbose() > 1) {
         printf("Boundary.minima_: %lf %lf %lf\n", 
                 boundary.minima_[0], boundary.minima_[1], boundary.minima_[2]);
         printf("Boundary.maxima_: %lf %lf %lf\n", 
                 boundary.maxima_[0], boundary.maxima_[1], boundary.maxima_[2]);
         printf("Boundary.L   : %lf %lf %lf\n", 
                 boundary.lengths_[0], boundary.lengths_[1], boundary.lengths_[2]);
      }
      #endif

      std::cout << std::endl;
      std::cout << "BravaisBasis(1)   " << boundary.bravaisBasisVector(1) << std::endl;
      std::cout << "ReciprocalBasis(1)" << boundary.reciprocalBasisVector(1) << std::endl;

   }

   void testStreamIO1() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/OrthorhombicBoundary", in);
      in >> boundary;
      in.close();

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testStreamIO2() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/TetragonalBoundary", in);
      in >> boundary;
      in.close();

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      TEST_ASSERT(eq(boundary.length(0), 2.0));
      TEST_ASSERT(eq(boundary.length(1), 3.0));
      TEST_ASSERT(eq(boundary.length(2), 3.0));

      Vector Lp = boundary.lengths();
      std::cout << "Lp = " << Lp << std::endl;
      TEST_ASSERT(eq(Lp[0], 2.0));
      TEST_ASSERT(eq(Lp[1], 3.0));
      TEST_ASSERT(eq(Lp[2], 3.0));

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testStreamIO3() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/CubicBoundary", in);
      in >> boundary;
      in.close();

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testSerialize() 
   {
      printMethod(TEST_FUNC);

      MemoryOArchive oar;
      MemoryIArchive iar;

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/TetragonalBoundary", in);
     
      in >> boundary;
      oar.allocate(2000);
      oar << boundary;
      iar = oar;

      OrthorhombicBoundary clone;
      iar >> clone;

      std::cout << std::endl;
      std::cout << clone << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < Dimension; i++) {
         TEST_ASSERT(eq(clone.minima_[i], 0.0));
         TEST_ASSERT(eq(clone.maxima_[i], clone.lengths_[i]));
         TEST_ASSERT(eq(clone.lengths_[i], boundary.lengths_[i]));
         TEST_ASSERT(eq(clone.volume(), boundary.volume()));
         TEST_ASSERT(clone.lengths_[i] > 1.0E-8);
      }

      #if 0
      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }
      #endif

   }

   void testShift()
   {
      printMethod(TEST_FUNC);

      Vector R, L;

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      boundary.setOrthorhombic(L);

      R[0] =  2.6;
      R[1] = -0.4;
      R[2] =  2.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 0.6));
      TEST_ASSERT(eq(R[1], 2.6));
      TEST_ASSERT(eq(R[2], 2.1));

      Vector Rg, Rn;
      boundary.transformCartToGen(R, Rg);
      boundary.transformGenToCart(Rg, Rn);
      for (int i = 0; i < Dimension; ++i) {
         TEST_ASSERT(Rg[i] >= 0.0);
         TEST_ASSERT(Rg[i] <  1.0);
         TEST_ASSERT(eq(R[i], Rn[i]));
      }

   };

   void testDistanceSq1()
   {
      printMethod(TEST_FUNC);

      Vector L;
      Vector R1;
      Vector R2;
      double dRSq;

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      boundary.setOrthorhombic(L);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  3.2;
      R2[0] =  1.2;
      R2[1] =  2.4;
      R2[2] =  0.1;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, 1.97));

   };

   void testDistanceSq2()
   {
      Vector    L;
      Vector    R1;
      Vector    R2;
      Vector    dR;
      IntVector shift; 
      double    dRSq1, dRSq2, dRSq3;

      printMethod(TEST_FUNC);

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      boundary.setOrthorhombic(L);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  3.2;
      R2[0] =  1.2;
      R2[1] =  2.4;
      R2[2] =  0.1;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);

      TEST_ASSERT(eq(dR[0], 0.4));
      TEST_ASSERT(eq(dR[1], 1.0));
      TEST_ASSERT(eq(dR[2], -0.9));
      TEST_ASSERT(eq(dRSq1, 1.97));
      TEST_ASSERT(eq(dRSq2, 1.97));
      TEST_ASSERT(eq(dRSq3, 1.97));

      R1[0] =  1.2;
      R1[1] =  2.4;
      R1[2] =  0.1;
      R2[0] =  1.6;
      R2[1] =  0.4;
      R2[2] =  3.2;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], - 0.4));
      TEST_ASSERT(eq(dR[1], - 1.0));
      TEST_ASSERT(eq(dR[2], 0.9));
      TEST_ASSERT(eq(dRSq1, 1.97));
      TEST_ASSERT(eq(dRSq2, 1.97));
      TEST_ASSERT(eq(dRSq3, 1.97));

   };

   void testTransforms()
   {
      printMethod(TEST_FUNC);

      Vector R, L;

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      boundary.setOrthorhombic(L);

      Random random;
      Vector Rg;
      Vector Rc;
      Vector Rg2;
      int i, j;
      random.setSeed(1793467182);
      for (i = 0; i < 1000; ++i) {
         for (j=0; j < Dimension; ++j) {
            Rg[j] = random.uniform(-0.5, 1.5);
         }
      }
      boundary.transformGenToCart(Rg, Rc);
      boundary.transformCartToGen(Rc, Rg2);
      
      for (int i = 0; i < Dimension; ++i) {
         TEST_ASSERT(Rg[i] >= 0.0);
         TEST_ASSERT(Rg[i] <  1.0);
         TEST_ASSERT(eq(Rg[i], Rg2[i]));
      }

   };

   void testInitializeCubic() 
   {
      printMethod(TEST_FUNC);

      Vector Lp;

      boundary.setCubic(2.0);

      // Assertions
      TEST_ASSERT(boundary.isValid());

      Lp = boundary.lengths();
      TEST_ASSERT(eq(Lp[0], 2.0));
      TEST_ASSERT(eq(Lp[1], 2.0));
      TEST_ASSERT(eq(Lp[2], 2.0));

      TEST_ASSERT(eq(boundary.length(0), 2.0));
      TEST_ASSERT(eq(boundary.length(1), 2.0));
      TEST_ASSERT(eq(boundary.length(2), 2.0));

   }

   void testStreamIO1Cubic() 
   {
      printMethod(TEST_FUNC);

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/CubicBoundary", in);
      in >> boundary;
      in.close();

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());

      Vector Lp = boundary.lengths();
      TEST_ASSERT(eq(Lp[0], 3.0));
      TEST_ASSERT(eq(Lp[1], 3.0));
      TEST_ASSERT(eq(Lp[2], 3.0));

      TEST_ASSERT(eq(boundary.length(0), 3.0));
      TEST_ASSERT(eq(boundary.length(1), 3.0));
      TEST_ASSERT(eq(boundary.length(2), 3.0));

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testShiftCubic()
   {
      printMethod(TEST_FUNC);

      Vector R;

      // Setup Boundary
      boundary.setCubic(2.0);

      R[0] =  1.6;
      R[1] = -0.6;
      R[2] =  2.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 1.6));
      TEST_ASSERT(eq(R[1], 1.4));
      TEST_ASSERT(eq(R[2], 0.1));

   }

   void testDistanceSq1Cubic()
   {
      printMethod(TEST_FUNC);

      Vector R1;
      Vector R2;
      double dRSq, result;

      // Setup Boundary
      boundary.setCubic(2.0);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  1.8;

      R2[0] =  1.2;
      R2[1] =  1.8;
      R2[2] =  0.1;

      result = 0.4*0.4 + 0.6*0.6 + 0.3*0.3;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, result));

   }

   void testDistanceSq2Cubic()
   {
      Vector R1;
      Vector R2;
      Vector dR; 
      double dRSq;

      printMethod(TEST_FUNC);

      // Setup Boundary
      boundary.setCubic(2.0);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  1.8;

      R2[0] =  1.2;
      R2[1] =  1.8;
      R2[2] =  0.1;
      dRSq = boundary.distanceSq(R1, R2, dR);

      TEST_ASSERT(eq(dR[0],  0.4));
      TEST_ASSERT(eq(dR[1],  0.6));
      TEST_ASSERT(eq(dR[2], -0.3));
      TEST_ASSERT(eq(dRSq, dR.square()));

      R1[0] =  1.2;
      R1[1] =  1.8;
      R1[2] =  0.1;

      R2[0] =  1.6;
      R2[1] =  0.4;
      R2[2] =  1.8;

      dRSq = boundary.distanceSq(R1, R2, dR);
      TEST_ASSERT(eq(dR[0], -0.4));
      TEST_ASSERT(eq(dR[1], -0.6));
      TEST_ASSERT(eq(dR[2],  0.3));
      TEST_ASSERT(eq(dRSq, dR.square()));

   }

   void testInitializeTetragonal() 
   {
      printMethod(TEST_FUNC);

      Vector Lp;

      boundary.setTetragonal(2.0, 3.0);

      // Assertions
      TEST_ASSERT(boundary.isValid());

      Lp = boundary.lengths();
      TEST_ASSERT(eq(Lp[0], 2.0));
      TEST_ASSERT(eq(Lp[1], 3.0));
      TEST_ASSERT(eq(Lp[2], 3.0));

      TEST_ASSERT(eq(boundary.length(0), 2.0));
      TEST_ASSERT(eq(boundary.length(1), 3.0));
      TEST_ASSERT(eq(boundary.length(2), 3.0));

   }

   void testStreamIOTetragonal() 
   {
      printMethod(TEST_FUNC);

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/TetragonalBoundary", in);
      in >> boundary;
      in.close();

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      TEST_ASSERT(boundary.isValid());

      TEST_ASSERT(eq(boundary.length(0), 2.0));
      TEST_ASSERT(eq(boundary.length(1), 3.0));
      TEST_ASSERT(eq(boundary.length(2), 3.0));

      Vector Lp = boundary.lengths();
      std::cout << "Lp = " << Lp << std::cout;
      TEST_ASSERT(eq(Lp[0], 2.0));
      TEST_ASSERT(eq(Lp[1], 3.0));
      TEST_ASSERT(eq(Lp[2], 3.0));

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testShiftTetragonal()
   {
      printMethod(TEST_FUNC);

      boundary.setTetragonal(2.0, 3.0);

      Vector R;
      R[0] =  1.6;
      R[1] = -0.6;
      R[2] =  3.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 1.6));
      TEST_ASSERT(eq(R[1], 2.4));
      TEST_ASSERT(eq(R[2], 0.1));

   };

   void testDistanceSq1Tetragonal()
   {
      printMethod(TEST_FUNC);

      Vector R1;
      Vector R2;
      double dRSq, result;

      // Setup Boundary
      boundary.setTetragonal(2.0, 3.0);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  2.8;

      R2[0] =  1.2;
      R2[1] =  2.8;
      R2[2] =  0.1;

      result = 0.4*0.4 + 0.6*0.6 + 0.3*0.3;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, result));

   };

   void testDistanceSq2Tetragonal()
   {
      Vector R1;
      Vector R2;
      Vector dR; 
      double dRSq;

      printMethod(TEST_FUNC);

      boundary.setTetragonal(2.0, 3.0);

      R1[0] =  1.6;
      R1[1] =  0.4;
      R1[2] =  2.8;

      R2[0] =  1.2;
      R2[1] =  2.8;
      R2[2] =  0.1;
      dRSq = boundary.distanceSq(R1, R2, dR);

      TEST_ASSERT(eq(dR[0],  0.4));
      TEST_ASSERT(eq(dR[1],  0.6));
      TEST_ASSERT(eq(dR[2], -0.3));
      TEST_ASSERT(eq(dRSq, dR.square()));

      R1[0] =  1.2;
      R1[1] =  2.8;
      R1[2] =  0.1;

      R2[0] =  1.6;
      R2[1] =  0.4;
      R2[2] =  2.8;

      dRSq = boundary.distanceSq(R1, R2, dR);
      TEST_ASSERT(eq(dR[0], -0.4));
      TEST_ASSERT(eq(dR[1], -0.6));
      TEST_ASSERT(eq(dR[2],  0.3));
      TEST_ASSERT(eq(dRSq, dR.square()));

   }

};

TEST_BEGIN(OrthorhombicBoundaryTest)
TEST_ADD(OrthorhombicBoundaryTest, testInitialize)
TEST_ADD(OrthorhombicBoundaryTest, testStreamIO1)
TEST_ADD(OrthorhombicBoundaryTest, testStreamIO2)
TEST_ADD(OrthorhombicBoundaryTest, testStreamIO3)
TEST_ADD(OrthorhombicBoundaryTest, testSerialize)
TEST_ADD(OrthorhombicBoundaryTest, testShift)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq1)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq2)

TEST_ADD(OrthorhombicBoundaryTest, testInitializeCubic)
TEST_ADD(OrthorhombicBoundaryTest, testStreamIO1Cubic)
TEST_ADD(OrthorhombicBoundaryTest, testShiftCubic)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq1Cubic)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq2Cubic)

TEST_ADD(OrthorhombicBoundaryTest, testInitializeTetragonal)
TEST_ADD(OrthorhombicBoundaryTest, testStreamIOTetragonal)
TEST_ADD(OrthorhombicBoundaryTest, testShiftTetragonal)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq1Tetragonal)
TEST_ADD(OrthorhombicBoundaryTest, testDistanceSq2Tetragonal)

TEST_END(OrthorhombicBoundaryTest)

#endif
