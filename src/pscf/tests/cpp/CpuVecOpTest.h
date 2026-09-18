#ifndef PSCF_CPU_VEC_OP_TEST_H
#define PSCF_CPU_VEC_OP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/VecOpCx.h>
#include <pscf/backend/cpp/complex.h>

#include <util/math/Constants.h>
#include <util/containers/DArray.h>

#include <complex>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CpuVecOpTest : public UnitTest
{

private:

   // Error tolerance for array equality
   constexpr static double tolerance_ = 1E-10;

   // Array size, large enough to require multiple blocks
   const static int n = 2048;

   // Input and output arrays, real and complex
   DArray<double> inReal, inReal2, outReal, outReal2;

   DArray<fftw_complex> inComplex, inComplex2, outComplex;

   // Input scalars, real and complex
   double scalarReal;
   fftw_complex scalarComplex;

   // Input and outputs using standard types, for comparison
   DArray<double> refOutReal, refOutReal2;
   DArray<fftw_complex> refOutComplex;

public:

   void setUp()
   {
      // Allocate arrays
      inReal.allocate(n);
      inReal2.allocate(n);
      outReal.allocate(n);
      outReal2.allocate(n);

      inComplex.allocate(n);
      inComplex2.allocate(n);
      outComplex.allocate(n);

      refOutReal.allocate(n);
      refOutReal2.allocate(n);
      refOutComplex.allocate(n);

      // Define "in" arrays with arbitrary data between -1 and 1
      double twoPi = 2.0 * Constants::Pi;
      double fourPi = 4.0 * Constants::Pi;
      for (int i = 0; i < n; i++) {
         double frac = (double)i / (double)n;

         inReal[i] = sin(fourPi * frac);
         inReal2[i] = cos(frac); // all values >0.5, for dividing
         inComplex[i][0] = cos(twoPi * frac);
         inComplex[i][1] = sin(twoPi * frac);
         inComplex2[i][0] = erf(frac - 0.5);
         inComplex2[i][1] = 1 - cos(frac);

         // note: the above two lines use copy-list-initialization
      }

      // Define input scalars with arbitrary values
      scalarReal = 0.633;
      scalarComplex[0] = -0.807;
      scalarComplex[1] = 0.0459;
   }

   void tearDown()
   {}

   // Test VecOp::eqV and VecOp::eqS
   void testEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test eqV (real) ~~~~~
      VecOp::eqV(outReal, inReal);
      checkEqualReal(outReal, inReal);

      VecOp::eqV(outComplex, inComplex);
      checkEqualComplex(outComplex, inComplex);

      // ~~~ Test eqS ~~~
      VecOp::eqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqS(outComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);
   }

   // Test VecOp::addVV and VecOp::addVS
   void testAdd()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addVV ~~~
      VecOp::addVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] + inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::addVV(outComplex, inComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         add(refOutComplex[i], inComplex[i], inComplex2[i]);
         //refOutComplex[i][0] = inComplex[i][0] + inComplex2[i][0];
         //refOutComplex[i][1] = inComplex[i][1] + inComplex2[i][1];
      }
      checkEqualComplex(outComplex, refOutComplex);

      VecOp::addVV(outComplex, inComplex, inReal);
      for (int i = 0; i < n; i++) {
         add(refOutComplex[i], inComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test addVS ~~~
      VecOp::addVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         add(refOutReal[i], inReal[i], scalarReal);
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::addVS(outComplex, inComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         add(refOutComplex[i], inComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      VecOp::addVS(outComplex, inComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         add(refOutComplex[i], inComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::subVV and VecOp::subVS
   void testSub()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subVV ~~~
      VecOp::subVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] - inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::subVV(outComplex, inComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         sub(refOutComplex[i], inComplex[i], inComplex2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      VecOp::subVV(outComplex, inComplex, inReal);
      for (int i = 0; i < n; i++) {
         sub(refOutComplex[i], inComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test subVS (real) ~~~
      VecOp::subVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] - scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test subVS (complex) ~~~
      VecOp::subVS(outComplex, inComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         sub(refOutComplex[i], inComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test subVS (mixed) ~~~
      VecOp::subVS(outComplex, inComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         sub(refOutComplex[i], inComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::mulVV and VecOp::mulVS
   void testMul()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulVV (real) ~~~
      VecOp::mulVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test mulVV (complex) ~~~
      VecOp::mulVV(outComplex, inComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         mul(refOutComplex[i], inComplex[i], inComplex2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test mulVV (mixed) ~~~
      VecOp::mulVV(outComplex, inComplex, inReal);
      for (int i = 0; i < n; i++) {
         mul(refOutComplex[i], inComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test mulVS (real) ~~~
      VecOp::mulVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test mulVS (complex) ~~~
      VecOp::mulVS(outComplex, inComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         mul(refOutComplex[i], inComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test mulVS (mixed) ~~~
      VecOp::mulVS(outComplex, inComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         mul(refOutComplex[i], inComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::divVV and VecOp::divVS
   void testDiv()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divVV ~~~
      VecOp::divVV(outReal, inReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] / inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::divVV(outComplex, inComplex, inReal2);
      for (int i = 0; i < n; i++) {
         div( refOutComplex[i], inComplex[i], inReal2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test divVS (real) ~~~
      VecOp::divVS(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] / scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test divVS (mixed) ~~~
      VecOp::divVS(outComplex, inComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         div(refOutComplex[i], inComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::addEqV and VecOp::addEqS
   void testAddEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addEqV (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::addEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] += inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test addEqV (complex) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::addEqV(outComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         addEq(refOutComplex[i], inComplex2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test addEqV (complex, real and imaginar inputs) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::addEqV(outComplex, inReal, inReal2);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         refOutComplex[i][0] += inReal[i];
         refOutComplex[i][1] += inReal2[i];
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test addEqV (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::addEqV(outComplex, inReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         addEq(refOutComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test addEqS (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::addEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] += scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test addEqS (complex) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::addEqS(outComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         addEq(refOutComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test addEqS (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::addEqS(outComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         addEq(refOutComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::subEqV and VecOp::subEqS
   void testSubEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test subEqV (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::subEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] -= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test subEqV (complex) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::subEqV(outComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         subEq(refOutComplex[i], inComplex2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test subEqV (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::subEqV(outComplex, inReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         subEq(refOutComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test subEqS (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::subEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] -= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test subEqS (complex) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::subEqS(outComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         subEq(refOutComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test subEqS (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::subEqS(outComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         subEq(refOutComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);
   }

   // Test VecOp::mulEqV and VecOp::mulEqS
   void testMulEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test mulEqV ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::mulEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] *= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::eqV(outComplex, inComplex);
      VecOp::mulEqV(outComplex, inComplex2);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         mulEq(refOutComplex[i], inComplex2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      VecOp::eqV(outComplex, inComplex);
      VecOp::mulEqV(outComplex, inReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         mulEq(refOutComplex[i], inReal[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test mulEqS (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::mulEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] *= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test mulEqS (complex) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::mulEqS(outComplex, scalarComplex);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         mulEq(refOutComplex[i], scalarComplex);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test mulEqS (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::mulEqS(outComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         mulEq(refOutComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::divEqV and VecOp::divEqS
   void testDivEq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test divEqV (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::divEqV(outReal, inReal2);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] /= inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test divEqV (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::divEqV(outComplex, inReal2);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         divEq(refOutComplex[i], inReal2[i]);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test divEqS (real) ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::divEqS(outReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] /= scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test divEqS (mixed) ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::divEqS(outComplex, scalarReal);
      for (int i = 0; i < n; i++) {
         assign(refOutComplex[i], inComplex[i]);
         divEq(refOutComplex[i], scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test VecOp::expV
   void testExp()
   {
      printMethod(TEST_FUNC);

      // ~~~ expV (real) ~~~~~~
      VecOp::expV(outReal, inReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = std::exp(inReal[i]);
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ expV (complex) ~~~~~~~
      VecOp::expV(outComplex, inComplex);
      std::complex<double> c, d;
      for (int i = 0; i < n; i++) {
         assign(c, inComplex[i]);
	 d = std::exp(c);
         assign(refOutComplex[i], d);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ expVc (real) ~~~~~~
      VecOp::expVc(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = std::exp(inReal[i]*scalarReal);
      }
      checkEqualReal(outReal, refOutReal);

   }

   // Test VecOp::sqV
   void testSq()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test using full arrays ~~~
      double vr;
      VecOp::sqV(outReal, inReal);
      for (int i = 0; i < n; i++) {
         vr = inReal[i];
         refOutReal[i] = vr * vr;
      }
      checkEqualReal(outReal, refOutReal);

      std::complex<double> c, d;
      VecOp::sqV(outComplex, inComplex);
      for (int i = 0; i < n; i++) {
         assign(c, inComplex[i]);
	 d = c*c;
         assign(refOutComplex[i], d);
      }
      checkEqualComplex(outComplex, refOutComplex);

   }

   // Test the other miscellaneous vector operations in Vec.h
   void testMisc()
   {
      printMethod(TEST_FUNC);

      // ~~~ Test addVcVc ~~~
      VecOp::addVcVc(outReal, inReal, scalarReal, inReal2, -2.0);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * scalarReal + (-2.0)*inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test addVcVcS ~~~
      VecOp::addVcVcS(outReal, inReal, scalarReal, inReal2, -2.0, 0.36);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * scalarReal + (-2.0)*inReal2[i] + 0.36;
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test addEqVc ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::addEqVc(outReal, inReal2, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal[i] += inReal2[i] * scalarReal;
      }
      checkEqualReal(outReal, refOutReal);

      #if 0
      // ~~~ Test subVVS ~~~
      VecOp::subVVS(outReal, inReal, inReal2, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] - inReal2[i] - scalarReal;
      }
      checkEqualReal(outReal, refOutReal);
      #endif

      #if 0
      // ~~~ Test divEqVc ~~~
      VecOp::eqV(outComplex, inComplex);
      VecOp::divEqVc(outComplex, inReal2, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutComplex[i] = inComplex[i];
         refOutComplex[i] /= (inReal2[i] * scalarReal);
      }
      checkEqualComplex(outComplex, refOutComplex);

      // ~~~ Test expVc ~~~
      VecOp::expVc(outReal, inReal, scalarReal);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = exp(inReal[i] * scalarReal);
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test eqVPair ~~~
      VecOp::eqVPair(outReal, outReal2, inReal);
      outReal2 = outReal2;
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i];
         refOutReal2[i] = inReal[i];
      }
      checkEqualReal(outReal, refOutReal);
      checkEqualReal(outReal2, refOutReal2);

      // ~~~ Test mulVVPair ~~~
      VecOp::mulVVPair(outReal, outReal2, inReal, inReal2, inReal2);
      outReal2 = outReal2;
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * inReal2[i];
         refOutReal2[i] = inReal2[i] * inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);
      checkEqualReal(outReal2, refOutReal2);

      // ~~~ Test mulEqVPair ~~~
      VecOp::eqV(outReal, inReal);
      VecOp::eqV(outReal2, inReal2);
      VecOp::mulEqVPair(outReal, outReal2, inReal2);
      outReal2 = outReal2;
      checkEqualReal(outReal, refOutReal);   // same ref. array as above
      checkEqualReal(outReal2, refOutReal2); // same ref. array as above

      // ~~~ Test addVMany ~~~
      DArray< DArray<double> const *> inVecs;
      inVecs.allocate(4);
      inVecs[0] = &inReal;
      inVecs[1] = &inReal2;
      inVecs[2] = &inReal;
      inVecs[3] = &inReal2;
      VecOp::addVMany(outReal, inVecs);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = 2 * (inReal[i] + inReal2[i]);
      }
      checkEqualReal(outReal, refOutReal);

      DArray< DArray<double> > inVecs2;
      inVecs2.allocate(4);
      inVecs2[0].associate(inReal, 0, inReal.capacity());
      inVecs2[1].associate(inReal2, 0, inReal2.capacity());
      inVecs2[2].associate(inReal, 0, inReal.capacity());
      inVecs2[3].associate(inReal2, 0, inReal2.capacity());
      VecOp::addVMany(outReal, inVecs2);
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test mulVMany ~~~
      VecOp::mulVMany(outReal, inVecs);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = inReal[i] * inReal[i] *
                         inReal2[i] * inReal2[i];
      }
      checkEqualReal(outReal, refOutReal);

      VecOp::mulVMany(outReal, inVecs2);
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test sqAbsV ~~~
      VecOp::sqAbsV(outReal, inComplex);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = norm(inComplex[i]);
      }
      checkEqualReal(outReal, refOutReal);

      // ~~~ Test sqSqAbsV ~~~
      VecOp::sqSqAbsV(outReal, inComplex);
      for (int i = 0; i < n; i++) {
         refOutReal[i] = std::pow(std::norm(inComplex[i]), 2.0);
      }
      checkEqualReal(outReal, refOutReal);
      #endif
   }

   void checkEqualReal(DArray<double>& a, DArray<double>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);

      for (int i = 0; i < n; i++) {
         TEST_ASSERT(std::abs(a[i] - b[i]) < tolerance_);
      }
   }

   void checkEqualComplex(DArray<fftw_complex>& a,
                          DArray<fftw_complex>& b)
   {
      int n = a.capacity();
      TEST_ASSERT(b.capacity() == n);
      TEST_ASSERT(n > 0);
      for (int i = 0; i < n; i++) {
         TEST_ASSERT(std::abs(a[i][0] - b[i][0]) < tolerance_);
         TEST_ASSERT(std::abs(a[i][1] - b[i][1]) < tolerance_);
      }
   }

};

TEST_BEGIN(CpuVecOpTest)
TEST_ADD(CpuVecOpTest, testEq)
TEST_ADD(CpuVecOpTest, testAdd)
TEST_ADD(CpuVecOpTest, testSub)
TEST_ADD(CpuVecOpTest, testMul)
TEST_ADD(CpuVecOpTest, testDiv)
TEST_ADD(CpuVecOpTest, testAddEq)
TEST_ADD(CpuVecOpTest, testSubEq)
TEST_ADD(CpuVecOpTest, testMulEq)
TEST_ADD(CpuVecOpTest, testDivEq)
TEST_ADD(CpuVecOpTest, testExp)
TEST_ADD(CpuVecOpTest, testSq)
TEST_ADD(CpuVecOpTest, testMisc)

TEST_END(CpuVecOpTest)

#endif
