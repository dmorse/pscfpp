#ifndef RPG_SOLVER_KERNEL_TEST_H
#define RPG_SOLVER_KERNEL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/solvers/Block.h>
#include <rpg/solvers/WaveList.h>
#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class SolverKernelTest : public UnitTest
{
   
private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   constexpr static float tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static double tolerance_ = 1E-10;
   #endif
   #endif

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testRealMulVConjVV()
   {
      printMethod(TEST_FUNC);

      int n = 100000;

      HostDArray<cudaReal> a_h(n), a_ref(n), d_h(n);
      HostDArray<cudaComplex> b_h(n), c_h(n);

      double twoPi = 2.0 * Constants::Pi;
      double fourPi = 4.0 * Constants::Pi;

      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         b_h[i].x = sin(fourPi * frac);
         b_h[i].y = cos(frac);
         c_h[i].x = cos(twoPi * frac);
         c_h[i].y = sin(twoPi * frac);
         d_h[i] = erf(frac - 0.5);

         // Calculate a on host
         a_ref[i] = ((b_h[i].x * c_h[i].x) + (b_h[i].y * c_h[i].y)) * d_h[i];
      }

      DeviceArray<cudaReal> a(n), d(n);
      DeviceArray<cudaComplex> b(n), c(n);
      b = b_h;
      c = c_h;
      d = d_h;

      // Perform calculation on device
      realMulVConjVV(a, b, c, d);
      a_h = a;

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a_ref[i] - a_h[i]) < tolerance_);
      }
   }

   void testMulVVPair()
   {
      printMethod(TEST_FUNC);

      int n = 100000;

      HostDArray<cudaReal> out1_h(n), out2_h(n), out1_ref(n), out2_ref(n), 
                           shared_h(n), in1_h(n), in2_h(n); 
      
      double fourPi = 4.0 * Constants::Pi;

      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         shared_h[i] = sin(fourPi * frac);
         in1_h[i] = cos(frac);
         in2_h[i] = erf(frac - 0.5);

         // Calculate out1 and out2 on host
         out1_ref[i] = shared_h[i] * in1_h[i];
         out2_ref[i] = shared_h[i] * in2_h[i];
      }

      DeviceArray<cudaReal> out1(n), out2(n), shared(n), in1(n), in2(n); 
      shared = shared_h;
      in1 = in1_h;
      in2 = in2_h;

      // Perform calculation on device
      mulVVPair(out1, out2, shared, in1, in2);
      out1_h = out1;
      out2_h = out2;

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(out1_h[i] - out1_ref[i]) < tolerance_);
         TEST_ASSERT(abs(out2_h[i] - out2_ref[i]) < tolerance_);
      }
   }

   void testMulEqVPair()
   {
      printMethod(TEST_FUNC);

      int n = 100000;

      HostDArray<cudaReal> out1_h(n), out2_h(n), out1_ref(n), out2_ref(n), 
                           shared_h(n); 
      
      double fourPi = 4.0 * Constants::Pi;

      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         shared_h[i] = sin(fourPi * frac);
         out1_h[i] = cos(frac);
         out2_h[i] = erf(frac - 0.5);

         // Perform calculation on host
         out1_ref[i] = shared_h[i] * out1_h[i];
         out2_ref[i] = shared_h[i] * out2_h[i];
      }

      DeviceArray<cudaReal> out1(n), out2(n), shared(n); 
      shared = shared_h;
      out1 = out1_h;
      out2 = out2_h;

      // Perform calculation on device
      mulEqVPair(out1, out2, shared);
      out1_h = out1;
      out2_h = out2;

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(out1_h[i] - out1_ref[i]) < tolerance_);
         TEST_ASSERT(abs(out2_h[i] - out2_ref[i]) < tolerance_);
      }
   }

   void testRichardsonEx()
   {
      printMethod(TEST_FUNC);

      int n = 100000;

      HostDArray<cudaReal> qNew_h(n), qNew_ref(n), 
                           qr_h(n), qr2_h(n), expW2_h(n);
      
      double fourPi = 4.0 * Constants::Pi;

      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         qr_h[i] = sin(fourPi * frac);
         qr2_h[i] = cos(frac);
         expW2_h[i] = erf(frac - 0.5);

         // Calculate qNew on host
         qNew_ref[i] = (4.0 * (qr2_h[i] * expW2_h[i]) - qr_h[i]) / 3.0;
      }

      DeviceArray<cudaReal> qNew(n), qr(n), qr2(n), expW2;
      qr = qr_h;
      qr2 = qr2_h;
      expW2 = expW2_h;

      // Perform calculation on device
      richardsonEx(qNew, qr, qr2, expW2);
      qNew_h = qNew;

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(qNew_h[i] - qNew_ref[i]) < tolerance_);
      }
   }

   void testAddEqMulVVc()
   {
      printMethod(TEST_FUNC);

      int n = 100000;

      HostDArray<cudaReal> a_h(n), a_ref(n), b_h(n), c_h(n);
      cudaReal d = -4.356;
      
      double fourPi = 4.0 * Constants::Pi;

      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         a_h[i] = sin(fourPi * frac);
         b_h[i] = cos(frac);
         c_h[i] = erf(frac - 0.5);

         // Calculate a on host
         a_ref[i] = a_h[i] + (b_h[i] * c_h[i] * d);
      }

      DeviceArray<cudaReal> a(n), b(n), c(n);
      a = a_h;
      b = b_h;
      c = c_h;

      // Perform calculation on device
      addEqMulVVc(a, b, c, d);
      a_h = a;

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a_h[i] - a_ref[i]) < tolerance_);
      }
   }
};

TEST_BEGIN(SolverKernelTest)
TEST_ADD(SolverKernelTest, testRealMulVConjVV)
TEST_ADD(SolverKernelTest, testMulVVPair)
TEST_ADD(SolverKernelTest, testMulEqVPair)
TEST_ADD(SolverKernelTest, testRichardsonEx)
TEST_ADD(SolverKernelTest, testAddEqMulVVc)
TEST_END(SolverKernelTest)

#endif