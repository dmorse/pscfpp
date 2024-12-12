#ifndef PSCF_CUDA_VEC_TEST_H
#define PSCF_CUDA_VEC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/Vec.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/DeviceDArray.h>
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/FieldComparison.h>
#include <util/math/Constants.h>
#include <complex>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CudaVecTest : public UnitTest
{

private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   typedef float numType;
   constexpr static numType tolerance = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   typedef double numType;
   constexpr static numType tolerance = 1E-10;
   #endif
   #endif

   // Array size, large enough to require multiple blocks
   const static int n = 2048; 

   // Input and output arrays, real and complex
   HostDArray<cudaReal> hInReal, hInReal2, hOutReal;
   DeviceDArray<cudaReal> dInReal, dInReal2, dOutReal;

   HostDArray<cudaComplex> hInComplex, hInComplex2, hOutComplex;
   DeviceDArray<cudaComplex> dInComplex, dInComplex2, dOutComplex;

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

   // Test Vec::eqV and Vec::eqS
   void testEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test eqV ~~~
      Vec::eqV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS ~~~
      Vec::eqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqV using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal2, n/4, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex2, n/4, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refInComplex[i+(n/4)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test eqS using partial arrays ~~~
      Vec::eqS(dOutReal, scalarReal*2, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutReal[i+(n/4)] = refScalarReal*2;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqS(dOutComplex, scalarComplex, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // edit reference array
         refOutComplex[i+(n/4)] = refScalarComplex*2;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::addVV and Vec::addVS
   void testAdd()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addVV ~~~
      Vec::addVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::addVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::addVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::addVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS ~~~
      Vec::addVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::addVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::addVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::addVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addVV using partial arrays ~~~
      Vec::addVV(dOutReal, dInReal, dInReal, n/2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = refInReal[i+(n/4)] + refInReal[i+(n/2)];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::addVV(dOutComplex, dInComplex, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex[i+(n/4)] + refInComplex[i+(n/2)];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::addVV(dOutComplex, dInReal2, dInComplex, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/4)] + refInComplex[i+(n/2)];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::addVV(dOutComplex, dInComplex, dInReal2, n/2, 0, n/2, n/4);
      hOutComplex = dOutComplex; 
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test addVS using partial arrays ~~~
      Vec::addVS(dOutReal, dInReal2, scalarReal, n/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i] + refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::addVS(dOutComplex, dInComplex2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::addVS(dOutComplex, dInReal2, scalarComplex, n/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i] + refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex, n/4, 0, n/4);

      Vec::addVS(dOutComplex, dInComplex2, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i] + refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::subVV and Vec::subVS
   void testSub()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subVV ~~~
      Vec::subVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::subVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::subVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS ~~~
      Vec::subVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::subVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVV using partial arrays ~~~
      Vec::subVV(dOutReal, dInReal2, dInReal, n/4, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/4)] - refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::subVV(dOutComplex, dInComplex2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVV(dOutComplex, dInReal2, dInComplex, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInReal2[i+(n/4)] - refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::subVV(dOutComplex, dInComplex, dInReal2, n/4, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex[i+(n/4)] - refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subVS using partial arrays ~~~
      Vec::subVS(dOutReal, dInReal2, scalarReal, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/2)] - refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::subVS(dOutComplex, dInComplex2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVS(dOutComplex, dInReal2, scalarComplex, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/2)] - refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::subVS(dOutComplex, dInComplex2, scalarReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/2)] - refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::mulVV and Vec::mulVS
   void testMul()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulVV ~~~
      Vec::mulVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::mulVV(dOutComplex, dInComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVV(dOutComplex, dInReal, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::mulVV(dOutComplex, dInComplex, dInReal);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS ~~~
      Vec::mulVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::mulVS(dOutComplex, dInComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVS(dOutComplex, dInReal, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInReal[i] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulVV using partial arrays ~~~
      Vec::mulVV(dOutReal, dInReal, dInReal, n/2, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/2)] = refInReal[i+(n/4)] * refInReal[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::mulVV(dOutComplex, dInComplex, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInComplex[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVV(dOutComplex, dInReal2, dInComplex, n/2, n/4, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/2)] = refInReal2[i+(n/4)] * refInComplex[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      Vec::mulVV(dOutComplex, dInComplex, dInReal2, n/2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test mulVS using partial arrays ~~~
      Vec::mulVS(dOutReal, dInReal2, scalarReal, 0, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::mulVS(dOutComplex, dInComplex2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVS(dOutComplex, dInReal2, scalarComplex, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInReal2[i+(n/4)] * refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::mulVS(dOutComplex, dInComplex2, scalarReal, 0, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] * refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::divVV and Vec::divVS
   void testDiv()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divVV ~~~
      Vec::divVV(dOutReal, dInReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::divVV(dOutComplex, dInComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS ~~~
      Vec::divVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::divVS(dOutComplex, dInComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divVV using partial arrays ~~~
      Vec::divVV(dOutReal, dInReal2, dInReal2, 0, n/4, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutReal[i] = refInReal2[i+(n/4)] / refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::divVV(dOutComplex, dInComplex2, dInReal2, 0, n/4, 0, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n*3/4; i++) { // get reference array
         refOutComplex[i] = refInComplex2[i+(n/4)] / refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex); 

      // ~~~ Test divVS using partial arrays ~~~
      Vec::divVS(dOutReal, dInReal2, scalarReal, n/4, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i+(n/4)] = refInReal2[i+(n/2)] / refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::divVS(dOutComplex, dInComplex2, scalarReal, n/4, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutComplex[i+(n/4)] = refInComplex2[i+(n/2)] / refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::expV
   void testExp()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test using full arrays ~~~
      Vec::expV(dOutReal, dInReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::expV(dOutComplex, dInComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = exp(refInComplex[i]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test using partial arrays ~~~
      Vec::expV(dOutReal, dInReal2, n/4, n/2, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutReal[i+(n/4)] = exp(refInReal2[i+(n/2)]);
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::expV(dOutComplex, dInComplex2, n/4, n/2, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n/4; i++) { // get reference array
         refOutComplex[i+(n/4)] = exp(refInComplex2[i+(n/2)]);
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::addEqV and Vec::addEqS
   void testAddEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addEqV ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::addEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::addEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] += refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqV using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::addEqV(dOutReal, dInReal2, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] += refInReal2[i+(n/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqV(dOutComplex, dInComplex2, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInComplex2[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqV(dOutComplex, dInReal, 0, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] += refInReal[i+(n/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test addEqS using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::addEqS(dOutReal, scalarReal, n/4, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/4) {
            refOutReal[i] += refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqS(dOutComplex, scalarComplex, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::addEqS(dOutComplex, scalarReal, n/4, n*3/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/4) {
            refOutComplex[i] += refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::subEqV and Vec::subEqS
   void testSubEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subEqV ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::subEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::subEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] -= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] -= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqV using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::subEqV(dOutReal, dInReal2, n/2, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n/2) {
            refOutReal[i] -= refInReal2[i-(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqV(dOutComplex, dInComplex2, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInComplex2[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqV(dOutComplex, dInReal, n/2, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n/2) {
            refOutComplex[i] -= refInReal[i-(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test subEqS using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::subEqS(dOutReal, scalarReal, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] -= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqS(dOutComplex, scalarComplex, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::subEqS(dOutComplex, scalarReal, 0, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] -= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::mulEqV and Vec::mulEqS
   void testMulEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulEqV ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::mulEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqV(dOutComplex, dInComplex2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInComplex2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqV(dOutComplex, dInReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refInReal[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::mulEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] *= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqS(dOutComplex, scalarComplex);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarComplex;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] *= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqV using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::mulEqV(dOutReal, dInReal2, 0, n/2, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n/2) {
            refOutReal[i] *= refInReal2[i+(n/2)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqV(dOutComplex, dInComplex2, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInComplex2[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqV(dOutComplex, dInReal, 0, n/2, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i < n/2) {
            refOutComplex[i] *= refInReal[i+(n/2)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test mulEqS using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::mulEqS(dOutReal, scalarReal, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutReal[i] *= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqS(dOutComplex, scalarComplex, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarComplex;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::mulEqS(dOutComplex, scalarReal, n/4, n/2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutComplex[i] *= refScalarReal;
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);
   }

   // Test Vec::divEqV and Vec::divEqS
   void testDivEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divEqV ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::divEqV(dOutReal, dInReal2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refInReal2[i];
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::divEqV(dOutComplex, dInReal2);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refInReal2[i];
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::divEqS(dOutReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] /= refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::divEqS(dOutComplex, scalarReal);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         refOutComplex[i] /= refScalarReal;
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqV using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::divEqV(dOutReal, dInReal2, n*3/4, 0, n/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i >= n*3/4) {
            refOutReal[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::divEqV(dOutComplex, dInReal2, n*3/4, 0, n/4);
      hOutComplex = dOutComplex;
      for (int i = 0; i < n; i++) { // get reference array
         refOutComplex[i] = refInComplex[i];
         if (i >= n*3/4) {
            refOutComplex[i] /= refInReal2[i-(n*3/4)];
         }
      }
      checkEqualComplex(hOutComplex, refOutComplex);

      // ~~~ Test divEqS using partial arrays ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::divEqS(dOutReal, scalarReal, 0, n*3/4);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if (i < n*3/4) {
            refOutReal[i] /= refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutComplex, dInComplex);
      Vec::divEqS(dOutComplex, scalarReal, 0, n*3/4);
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

      // ~~~ Test addEqMulVS ~~~
      Vec::eqV(dOutReal, dInReal);
      Vec::addEqMulVS(dOutReal, dInReal2, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         refOutReal[i] += refInReal2[i] * refScalarReal;
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::eqV(dOutReal, dInReal);
      Vec::addEqMulVS(dOutReal, dInReal2, scalarReal, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = refInReal[i];
         if ((i < n*3/4) && (i >= n/4)) {
            refOutReal[i] += refInReal2[i-(n/4)] * refScalarReal;
         }
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test sqNormV ~~~
      Vec::sqNormV(dOutReal, dInComplex);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = norm(refInComplex[i]);
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::sqNormV(dOutReal, dInComplex2, n/4, 0, n/2);
      hOutReal = dOutReal;
      for (int i = n/4; i < n*3/4; i++) { // get reference array
         refOutReal[i] = norm(refInComplex2[i-(n/4)]);
      }
      checkEqualReal(hOutReal, refOutReal);

      // ~~~ Test expMulVS ~~~
      Vec::expMulVS(dOutReal, dInReal, scalarReal);
      hOutReal = dOutReal;
      for (int i = 0; i < n; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i] * refScalarReal);
      }
      checkEqualReal(hOutReal, refOutReal);

      Vec::expMulVS(dOutReal, dInReal2, scalarReal, 0, n/4, n/2);
      hOutReal = dOutReal;
      for (int i = 0; i < n/2; i++) { // get reference array
         refOutReal[i] = exp(refInReal[i+(n/4)] * refScalarReal);
      }
      checkEqualReal(hOutReal, refOutReal);
   }

   void checkEqualReal(HostDArray<cudaReal>& a, DArray<numType>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i] - b[i]) < tolerance);
      }
   }

   void checkEqualComplex(HostDArray<cudaComplex>& a, 
                          DArray<std::complex<numType> >& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(abs(a[i].x - b[i].real()) < tolerance);
         TEST_ASSERT(abs(a[i].y - b[i].imag()) < tolerance);
      }
   }

};

TEST_BEGIN(CudaVecTest)
TEST_ADD(CudaVecTest, testEq)
TEST_ADD(CudaVecTest, testAdd)
TEST_ADD(CudaVecTest, testSub)
TEST_ADD(CudaVecTest, testMul)
TEST_ADD(CudaVecTest, testDiv)
TEST_ADD(CudaVecTest, testExp)
TEST_ADD(CudaVecTest, testAddEq)
TEST_ADD(CudaVecTest, testSubEq)
TEST_ADD(CudaVecTest, testMulEq)
TEST_ADD(CudaVecTest, testDivEq)
TEST_ADD(CudaVecTest, testMisc)
TEST_END(CudaVecTest)

#endif
