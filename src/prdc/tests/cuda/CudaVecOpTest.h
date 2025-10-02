#ifndef PRDC_CUDA_VEC_OP_TEST_H
#define PRDC_CUDA_VEC_OP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/resources.h>

#include <pscf/math/FieldComparison.h>
#include <util/math/Constants.h>
#include <complex>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

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
   HostDArray<Cuda::cudaReal> hInReal, hInReal2, hOutReal, hOutReal2;
   DeviceArray<Cuda::cudaReal> dInReal, dInReal2, dOutReal, dOutReal2;

   HostDArray<Cuda::cudaComplex> hInComplex, hInComplex2, hOutComplex;
   DeviceArray<Cuda::cudaComplex> dInComplex, dInComplex2, dOutComplex;

   // Input scalars, real and complex
   Cuda::cudaReal scalarReal;
   Cuda::cudaComplex scalarComplex;

   // Input and outputs using standard types, for comparison
   DArray<numType> refInReal, refInReal2, refOutReal, refOutReal2;
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
      hOutReal2.allocate(n);
      dInReal.allocate(n);
      dInReal2.allocate(n);
      dOutReal.allocate(n);
      dOutReal2.allocate(n);

      hInComplex.allocate(n);
      hInComplex2.allocate(n);
      hOutComplex.allocate(n);
      dInComplex.allocate(n);
      dInComplex2.allocate(n);
      dOutComplex.allocate(n);

      refInReal.allocate(n);
      refInReal2.allocate(n);
      refOutReal.allocate(n);
      refOutReal2.allocate(n);
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

   // Test Cuda::VecOp::eqV and Cuda::VecOp::eqS
   void testEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test eqV ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS ~~~
      Cuda::VecOp::eqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqV using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal2, n/4, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex2, n/4, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS using partial arrays ~~~
      Cuda::VecOp::eqS(dOutReal, scalarReal*2.0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refScalarReal*2.0;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::cudaComplex scalarComplex2;
      scalarComplex2.x = scalarComplex.x * 2.0;
      scalarComplex2.y = scalarComplex.y * 2.0;
      Cuda::VecOp::eqS(dOutComplex, scalarComplex2, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refScalarComplex * 2.0;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::addVV and Cuda::VecOp::addVS
   void testAdd()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addVV ~~~
      Cuda::VecOp::addVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::addVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::addVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS ~~~
      Cuda::VecOp::addVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::addVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addVV using partial arrays ~~~
      Cuda::VecOp::addVV(dOutReal, dInReal, dInReal, n/2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/2)] = refInReal[i] + refInReal[i+(n/4)];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::addVV(dOutComplex, dInComplex, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInComplex[i] + refInComplex[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVV(dOutComplex, dInReal2, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInReal2[i] + refInComplex[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::addVV(dOutComplex, dInComplex, dInReal2, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS using partial arrays ~~~
      Cuda::VecOp::addVS(dOutReal, dInReal2, scalarReal, n/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::addVS(dOutComplex, dInComplex2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVS(dOutComplex, dInReal2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::addVS(dOutComplex, dInComplex2, scalarReal, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::subVV and Cuda::VecOp::subVS
   void testSub()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subVV ~~~
      Cuda::VecOp::subVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::subVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::subVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS ~~~
      Cuda::VecOp::subVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::subVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVV using partial arrays ~~~
      Cuda::VecOp::subVV(dOutReal, dInReal2, dInReal, n/4, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)] - refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::subVV(dOutComplex, dInComplex2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVV(dOutComplex, dInReal2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::subVV(dOutComplex, dInComplex, dInReal2, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex[i+(n/4)] - refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS using partial arrays ~~~
      Cuda::VecOp::subVS(dOutReal, dInReal2, scalarReal, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/2)] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::subVS(dOutComplex, dInComplex2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVS(dOutComplex, dInReal2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::subVS(dOutComplex, dInComplex2, scalarReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::mulVV and Cuda::VecOp::mulVS
   void testMul()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulVV ~~~
      Cuda::VecOp::mulVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::mulVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::mulVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS ~~~
      Cuda::VecOp::mulVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::mulVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulVV using partial arrays ~~~
      Cuda::VecOp::mulVV(dOutReal, dInReal, dInReal, n/2, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/2)] = refInReal[i+(n/4)] * refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::mulVV(dOutComplex, dInComplex, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInComplex[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVV(dOutComplex, dInReal2, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInReal2[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Cuda::VecOp::mulVV(dOutComplex, dInComplex, dInReal2, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS using partial arrays ~~~
      Cuda::VecOp::mulVS(dOutReal, dInReal2, scalarReal, 0, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/4)] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::mulVS(dOutComplex, dInComplex2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVS(dOutComplex, dInReal2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::mulVS(dOutComplex, dInComplex2, scalarReal, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::divVV and Cuda::VecOp::divVS
   void testDiv()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divVV ~~~
      Cuda::VecOp::divVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::divVV(dOutComplex, dInComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS ~~~
      Cuda::VecOp::divVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::divVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divVV using partial arrays ~~~
      Cuda::VecOp::divVV(dOutReal, dInReal2, dInReal2, 0, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/4)] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::divVV(dOutComplex, dInComplex2, dInReal2, 0, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS using partial arrays ~~~
      Cuda::VecOp::divVS(dOutReal, dInReal2, scalarReal, n/4, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/2)] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::divVS(dOutComplex, dInComplex2, scalarReal, n/4, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/2)] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::expV
   void testExp()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test using full arrays ~~~
      Cuda::VecOp::expV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::expV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = exp(refInComplex[i]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test using partial arrays ~~~
      Cuda::VecOp::expV(dOutReal, dInReal2, n/4, n/2, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = exp(refInReal2[i+(n/2)]);
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::expV(dOutComplex, dInComplex2, n/4, n/2, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = exp(refInComplex2[i+(n/2)]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::addEqV and Cuda::VecOp::addEqS
   void testAddEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addEqV ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::addEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::addEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqV using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::addEqV(dOutReal, dInReal2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] += refInReal2[i+(n/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqV(dOutComplex, dInComplex2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInComplex2[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqV(dOutComplex, dInReal, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInReal[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::addEqS(dOutReal, scalarReal, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/4) {
            refOutReal[i] += refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqS(dOutComplex, scalarComplex, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::addEqS(dOutComplex, scalarReal, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::subEqV and Cuda::VecOp::subEqS
   void testSubEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subEqV ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::subEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::subEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqV using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::subEqV(dOutReal, dInReal2, n/2, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/2) {
            refOutReal[i] -= refInReal2[i-(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqV(dOutComplex, dInComplex2, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInComplex2[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqV(dOutComplex, dInReal, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInReal[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::subEqS(dOutReal, scalarReal, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] -= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqS(dOutComplex, scalarComplex, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::subEqS(dOutComplex, scalarReal, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::mulEqV and Cuda::VecOp::mulEqS
   void testMulEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulEqV ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::mulEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::mulEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqV using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::mulEqV(dOutReal, dInReal2, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] *= refInReal2[i+(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqV(dOutComplex, dInComplex2, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInComplex2[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqV(dOutComplex, dInReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInReal[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::mulEqS(dOutReal, scalarReal, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutReal[i] *= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqS(dOutComplex, scalarComplex, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::mulEqS(dOutComplex, scalarReal, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Cuda::VecOp::divEqV and Cuda::VecOp::divEqS
   void testDivEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divEqV ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::divEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::divEqV(dOutComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::divEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::divEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqV using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::divEqV(dOutReal, dInReal2, n*3/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n*3/4) {
            refOutReal[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::divEqV(dOutComplex, dInReal2, n*3/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n*3/4) {
            refOutComplex[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS using partial arrays ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::divEqS(dOutReal, scalarReal, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n*3/4) {
            refOutReal[i] /= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::divEqS(dOutComplex, scalarReal, 0, n*3/4);
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

      // ~~~ Test addVcVc ~~~
      Cuda::VecOp::addVcVc(dOutReal, dInReal, scalarReal, dInReal2, -1.0);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refScalarReal - refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test addVcVcVc ~~~
      Cuda::VecOp::addVcVcVc(dOutReal, dInReal, scalarReal, dInReal2, -1.0, dInReal, 1.0);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * (refScalarReal + 1) - refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test addEqVc ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::addEqVc(dOutReal, dInReal2, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test subVVS ~~~
      Cuda::VecOp::subVVS(dOutReal, dInReal, dInReal2, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refInReal2[i] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test divEqVc ~~~
      Cuda::VecOp::eqV(dOutComplex, dInComplex);
      Cuda::VecOp::divEqVc(dOutComplex, dInReal2, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= (refInReal2[i] * refScalarReal);
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test expVc ~~~
      Cuda::VecOp::expVc(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i] * refScalarReal);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test eqVPair ~~~
      Cuda::VecOp::eqVPair(dOutReal, dOutReal2, dInReal);
      hOutReal = dOutReal;
      hOutReal2 = dOutReal2;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal2[i] = refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);
      checkEqualReal(hOutReal2, refOutReal2);

      // ~~~ Test mulVVPair ~~~
      Cuda::VecOp::mulVVPair(dOutReal, dOutReal2, dInReal, dInReal2, dInReal2);
      hOutReal = dOutReal;
      hOutReal2 = dOutReal2;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal2[i];
         refOutReal2[i] = refInReal2[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);
      checkEqualReal(hOutReal2, refOutReal2);

      // ~~~ Test mulEqVPair ~~~
      Cuda::VecOp::eqV(dOutReal, dInReal);
      Cuda::VecOp::eqV(dOutReal2, dInReal2);
      Cuda::VecOp::mulEqVPair(dOutReal, dOutReal2, dInReal2);
      hOutReal = dOutReal;
      hOutReal2 = dOutReal2;
      checkEqualReal(hOutReal, refOutReal);   // same ref. array as above
      checkEqualReal(hOutReal2, refOutReal2); // same ref. array as above

      // ~~~ Test addVMany ~~~
      DArray<DeviceArray<Cuda::cudaReal> const *> inVecs;
      inVecs.allocate(4);
      inVecs[0] = &dInReal;
      inVecs[1] = &dInReal2;
      inVecs[2] = &dInReal;
      inVecs[3] = &dInReal2;
      Cuda::VecOp::addVMany(dOutReal, inVecs);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = 2 * (refInReal[i] + refInReal2[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      DArray<DeviceArray<Cuda::cudaReal> > inVecs2;
      inVecs2.allocate(4);
      inVecs2[0].associate(dInReal, 0, dInReal.capacity());
      inVecs2[1].associate(dInReal2, 0, dInReal2.capacity());
      inVecs2[2].associate(dInReal, 0, dInReal.capacity());
      inVecs2[3].associate(dInReal2, 0, dInReal2.capacity());
      Cuda::VecOp::addVMany(dOutReal, inVecs2);
      hOutReal = dOutReal;
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test mulVMany ~~~
      Cuda::VecOp::mulVMany(dOutReal, inVecs);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal[i] * 
                         refInReal2[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Cuda::VecOp::mulVMany(dOutReal, inVecs2);
      hOutReal = dOutReal;
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test sqNormV ~~~
      Cuda::VecOp::sqNormV(dOutReal, dInComplex);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = norm(refInComplex[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test sqSqNormV ~~~
      Cuda::VecOp::sqSqNormV(dOutReal, dInComplex);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = std::pow(std::norm(refInComplex[i]), 2.0);
      }
      checkEqualReal(hOutReal, refOutReal);
   }

   void checkEqualReal(HostDArray<Cuda::cudaReal>& a, DArray<numType>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i] - b[i]) < tolerance_);
      }
   }

   void checkEqualComplex(HostDArray<Cuda::cudaComplex>& a, 
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
