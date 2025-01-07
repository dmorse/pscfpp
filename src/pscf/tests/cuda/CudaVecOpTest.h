#ifndef PSCF_CUDA_VEC_OP_TEST_H
#define PSCF_CUDA_VEC_OP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/FieldComparison.h>
#include <util/math/Constants.h>
#include <complex>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CudaVecOpTest : public UnitTest
{

private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   typedef float numType;
   constexpr static numType tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   typedef double numType;
   constexpr static numType tolerance_ = 1E-10;
   #endif
   #endif

   // Array size, large enough to require multiple blocks
   const static int n = 2048; 

   // Input and output arrays, real and complex
   HostDArray<cudaReal> hInReal, hInReal2, hOutReal;
   DeviceArray<cudaReal> dInReal, dInReal2, dOutReal;

   HostDArray<cudaComplex> hInComplex, hInComplex2, hOutComplex;
   DeviceArray<cudaComplex> dInComplex, dInComplex2, dOutComplex;

   // Input scalars, real and complex
   cudaReal scalarReal;
   cudaComplex scalarComplex;

   // Input and outputs using standard types, for comparison
   DArray<numType> refInReal, refInReal2, refOutReal;
   DArray<std::complex<numType> > refInComplex, refInComplex2, 
                                  refOutComplex;
   numType refScalarReal;
   std::complex<numType> refScalarComplex;

public:

   void setUp()
   {
      // Allocate arrays
      hInReal.allocate(n);
      hInReal2.allocate(n);
      hOutReal.allocate(n);
      dInReal.allocate(n);
      dInReal2.allocate(n);
      dOutReal.allocate(n);

      hInComplex.allocate(n);
      hInComplex2.allocate(n);
      hOutComplex.allocate(n);
      dInComplex.allocate(n);
      dInComplex2.allocate(n);
      dOutComplex.allocate(n);

      refInReal.allocate(n);
      refInReal2.allocate(n);
      refOutReal.allocate(n);
      refInComplex.allocate(n);
      refInComplex2.allocate(n);
      refOutComplex.allocate(n);

      // Define "in" arrays with arbitrary data between -1 and 1
      double twoPi = 2.0 * Constants::Pi;
      double fourPi = 4.0 * Constants::Pi;
      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         hInReal[i] = sin(fourPi * frac);
         hInReal2[i] = cos(frac); // all values >0.5, for dividing
         hInComplex[i].x = cos(twoPi * frac);
         hInComplex[i].y = sin(twoPi * frac);
         hInComplex2[i].x = erf(frac - 0.5);
         hInComplex2[i].y = 1 - cos(frac);

         refInReal[i] = hInReal[i];
         refInReal2[i] = hInReal2[i];
         refInComplex[i] = {hInComplex[i].x, (hInComplex[i].y)};
         refInComplex2[i] = {hInComplex2[i].x, (hInComplex2[i].y)};
         // note: the above two lines use copy-list-initialization
      }

      // Copy from host to device
      dInReal = hInReal; 
      dInReal2 = hInReal2;
      dInComplex = hInComplex;;
      dInComplex2 = hInComplex2;

      // Define input scalars with arbitrary values
      scalarReal = 0.633; 
      scalarComplex.x = -0.807;
      scalarComplex.y = 0.0459;
      refScalarReal = scalarReal;
      refScalarComplex = {scalarComplex.x, scalarComplex.y};
   }

   void tearDown()
   {}

   // Test VecOp::eqV and VecOp::eqS
   void testEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test eqV ~~~
      VecOp::eqV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS ~~~
      VecOp::eqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqV using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal2, n/4, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex2, n/4, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS using partial arrays ~~~
      VecOp::eqS(dOutReal, scalarReal*2.0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refScalarReal*2.0;
      }
      checkEqualReal(hOutReal, refOutReal);

      cudaComplex scalarComplex2;
      scalarComplex2.x = scalarComplex.x * 2.0;
      scalarComplex2.y = scalarComplex.y * 2.0;
      VecOp::eqS(dOutComplex, scalarComplex2, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refScalarComplex * 2.0;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::addVV and VecOp::addVS
   void testAdd()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addVV ~~~
      VecOp::addVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::addVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::addVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS ~~~
      VecOp::addVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::addVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addVV using partial arrays ~~~
      VecOp::addVV(dOutReal, dInReal, dInReal, n/2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/2)] = refInReal[i] + refInReal[i+(n/4)];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::addVV(dOutComplex, dInComplex, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInComplex[i] + refInComplex[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVV(dOutComplex, dInReal2, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInReal2[i] + refInComplex[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::addVV(dOutComplex, dInComplex, dInReal2, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS using partial arrays ~~~
      VecOp::addVS(dOutReal, dInReal2, scalarReal, n/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::addVS(dOutComplex, dInComplex2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVS(dOutComplex, dInReal2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::addVS(dOutComplex, dInComplex2, scalarReal, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::subVV and VecOp::subVS
   void testSub()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subVV ~~~
      VecOp::subVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::subVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::subVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS ~~~
      VecOp::subVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::subVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVV using partial arrays ~~~
      VecOp::subVV(dOutReal, dInReal2, dInReal, n/4, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)] - refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::subVV(dOutComplex, dInComplex2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVV(dOutComplex, dInReal2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::subVV(dOutComplex, dInComplex, dInReal2, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex[i+(n/4)] - refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS using partial arrays ~~~
      VecOp::subVS(dOutReal, dInReal2, scalarReal, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/2)] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::subVS(dOutComplex, dInComplex2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVS(dOutComplex, dInReal2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::subVS(dOutComplex, dInComplex2, scalarReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::mulVV and VecOp::mulVS
   void testMul()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulVV ~~~
      VecOp::mulVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::mulVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::mulVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS ~~~
      VecOp::mulVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::mulVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulVV using partial arrays ~~~
      VecOp::mulVV(dOutReal, dInReal, dInReal, n/2, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/2)] = refInReal[i+(n/4)] * refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::mulVV(dOutComplex, dInComplex, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInComplex[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVV(dOutComplex, dInReal2, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInReal2[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      VecOp::mulVV(dOutComplex, dInComplex, dInReal2, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS using partial arrays ~~~
      VecOp::mulVS(dOutReal, dInReal2, scalarReal, 0, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/4)] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::mulVS(dOutComplex, dInComplex2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVS(dOutComplex, dInReal2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::mulVS(dOutComplex, dInComplex2, scalarReal, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::divVV and VecOp::divVS
   void testDiv()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divVV ~~~
      VecOp::divVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::divVV(dOutComplex, dInComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS ~~~
      VecOp::divVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::divVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divVV using partial arrays ~~~
      VecOp::divVV(dOutReal, dInReal2, dInReal2, 0, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/4)] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::divVV(dOutComplex, dInComplex2, dInReal2, 0, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS using partial arrays ~~~
      VecOp::divVS(dOutReal, dInReal2, scalarReal, n/4, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/2)] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::divVS(dOutComplex, dInComplex2, scalarReal, n/4, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/2)] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::expV
   void testExp()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test using full arrays ~~~
      VecOp::expV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::expV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = exp(refInComplex[i]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test using partial arrays ~~~
      VecOp::expV(dOutReal, dInReal2, n/4, n/2, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = exp(refInReal2[i+(n/2)]);
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::expV(dOutComplex, dInComplex2, n/4, n/2, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = exp(refInComplex2[i+(n/2)]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::addEqV and VecOp::addEqS
   void testAddEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addEqV ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqV using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqV(dOutReal, dInReal2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] += refInReal2[i+(n/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqV(dOutComplex, dInComplex2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInComplex2[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqV(dOutComplex, dInReal, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInReal[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqS(dOutReal, scalarReal, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/4) {
            refOutReal[i] += refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqS(dOutComplex, scalarComplex, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::addEqS(dOutComplex, scalarReal, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::subEqV and VecOp::subEqS
   void testSubEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subEqV ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::subEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::subEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqV using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::subEqV(dOutReal, dInReal2, n/2, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/2) {
            refOutReal[i] -= refInReal2[i-(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqV(dOutComplex, dInComplex2, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInComplex2[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqV(dOutComplex, dInReal, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInReal[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::subEqS(dOutReal, scalarReal, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] -= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqS(dOutComplex, scalarComplex, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::subEqS(dOutComplex, scalarReal, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::mulEqV and VecOp::mulEqS
   void testMulEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulEqV ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::mulEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::mulEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqV using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::mulEqV(dOutReal, dInReal2, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] *= refInReal2[i+(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqV(dOutComplex, dInComplex2, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInComplex2[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqV(dOutComplex, dInReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInReal[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::mulEqS(dOutReal, scalarReal, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutReal[i] *= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqS(dOutComplex, scalarComplex, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::mulEqS(dOutComplex, scalarReal, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test VecOp::divEqV and VecOp::divEqS
   void testDivEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divEqV ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::divEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::divEqV(dOutComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::divEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::divEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqV using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::divEqV(dOutReal, dInReal2, n*3/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n*3/4) {
            refOutReal[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::divEqV(dOutComplex, dInReal2, n*3/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n*3/4) {
            refOutComplex[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS using partial arrays ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::divEqS(dOutReal, scalarReal, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n*3/4) {
            refOutReal[i] /= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutComplex, dInComplex);
      VecOp::divEqS(dOutComplex, scalarReal, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n*3/4) {
            refOutComplex[i] /= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test the other miscellaneous vector operations in Vec.h
   void testMisc()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addEqVc ~~~
      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqVc(dOutReal, dInReal2, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::eqV(dOutReal, dInReal);
      VecOp::addEqVc(dOutReal, dInReal2, scalarReal, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutReal[i] += refInReal2[i-(n/4)] * refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test sqNormV ~~~
      VecOp::sqNormV(dOutReal, dInComplex);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = norm(refInComplex[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::sqNormV(dOutReal, dInComplex2, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = n/4; i < n*3/4; i++) { // get reference array
         refOutReal[i] = norm(refInComplex2[i-(n/4)]);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test expVc ~~~
      VecOp::expVc(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i] * refScalarReal);
      }
      checkEqualReal(hOutReal, refOutReal);

      VecOp::expVc(dOutReal, dInReal2, scalarReal, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = exp(refInReal2[i+(n/4)] * refScalarReal);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test addVMany ~~~
      DArray<DeviceArray<cudaReal> const *> inVecs;
      inVecs.allocate(4);
      inVecs[0] = &dInReal;
      inVecs[1] = &dInReal2;
      inVecs[2] = &dInReal;
      inVecs[3] = &dInReal2;
      VecOp::addVMany(dOutReal, inVecs);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = 2 * (refInReal[i] + refInReal2[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test mulVMany ~~~
      VecOp::mulVMany(dOutReal, inVecs);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal[i] * 
                         refInReal2[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);
   }

   void checkEqualReal(HostDArray<cudaReal>& a, DArray<numType>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i] - b[i]) < tolerance_);
      }
   }

   void checkEqualComplex(HostDArray<cudaComplex>& a, 
                          DArray<std::complex<numType> >& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i].x - b[i].real()) < tolerance_);
         TEST_ASSERT(abs(a[i].y - b[i].imag()) < tolerance_);
      }
   }

};

TEST_BEGIN(CudaVecOpTest)
TEST_ADD(CudaVecOpTest, testEq)
TEST_ADD(CudaVecOpTest, testAdd)
TEST_ADD(CudaVecOpTest, testSub)
TEST_ADD(CudaVecOpTest, testMul)
TEST_ADD(CudaVecOpTest, testDiv)
TEST_ADD(CudaVecOpTest, testExp)
TEST_ADD(CudaVecOpTest, testAddEq)
TEST_ADD(CudaVecOpTest, testSubEq)
TEST_ADD(CudaVecOpTest, testMulEq)
TEST_ADD(CudaVecOpTest, testDivEq)
TEST_ADD(CudaVecOpTest, testMisc)
TEST_END(CudaVecOpTest)

#endif
