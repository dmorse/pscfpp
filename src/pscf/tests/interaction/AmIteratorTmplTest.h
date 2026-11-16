#ifndef PSCF_AM_ITERATOR_TMPL_TEST_H
#define PSCF_AM_ITERATOR_TMPL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/iterator/AmIteratorTmpl.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

#include <pscf/iterator/AmIteratorTmpl.tpp>

#include <gsl/gsl_errno.h>

#include <vector>
#include <cmath>

using namespace Util;
using namespace Pscf;

namespace {

   // Number of calls of the GSL error handler
   int amIteratorTmplTestGslErrors = 0;

   // GSL error handler that counts errors instead of aborting
   void amIteratorTmplTestGslHandler(const char*, const char*, int, int)
   {  ++amIteratorTmplTestGslErrors; }

}

/*
* Fixed point problem for AmIteratorTmpl, with state and residual
* vectors of 4 elements. The first evaluations return prescribed
* residuals; later evaluations return the residual R = X0 - X, for
* which the solution is X = X0.
*/
class AmTestIterator
 : public AmIteratorTmpl< ParamComposite, DArray<double> >
{

public:

   AmTestIterator()
   {
      epsilon_ = 1.0E-10;
      maxItr_ = 100;
      maxHist_ = 10;
      verbose_ = 0;
      errorType_ = "rms";
      UTIL_CHECK(isValidErrorType());
      x_.allocate(4);
      r_.allocate(4);
      for (int i = 0; i < 4; ++i) {
         x_[i] = 0.0;
         r_[i] = 0.0;
      }
   }

   // Prescribed residuals of the first evaluations
   std::vector< std::vector<double> > residuals;

   // Largest absolute difference between the state and the solution
   double error() const
   {
      double error = 0.0;
      for (int i = 0; i < 4; ++i) {
         double diff = std::abs(x_[i] - double(i + 1));
         if (!(diff <= error)) error = diff; // also if diff is NaN
      }
      return error;
   }

private:

   DArray<double> x_;
   DArray<double> r_;
   int nEvaluate_ = 0;

   int nElements() override
   {  return 4; }

   bool hasInitialGuess() override
   {  return true; }

   void getCurrent(DArray<double>& curr) override
   {
      for (int i = 0; i < 4; ++i) curr[i] = x_[i];
   }

   void evaluate() override
   {
      for (int i = 0; i < 4; ++i) {
         if (nEvaluate_ < (int)residuals.size()) {
            r_[i] = residuals[nEvaluate_][i];
         } else {
            r_[i] = double(i + 1) - x_[i];
         }
      }
      ++nEvaluate_;
   }

   void getResidual(DArray<double>& resid) override
   {
      for (int i = 0; i < 4; ++i) resid[i] = r_[i];
   }

   void update(DArray<double>& newGuess) override
   {
      for (int i = 0; i < 4; ++i) x_[i] = newGuess[i];
   }

   void outputToLog() override
   {}

};

class AmIteratorTmplTest : public UnitTest
{

public:

   void setUp()
   {
      amIteratorTmplTestGslErrors = 0;
      oldHandler_ = gsl_set_error_handler(&amIteratorTmplTestGslHandler);
   }

   void tearDown()
   {  gsl_set_error_handler(oldHandler_); }

   /*
   * Iterations with an exactly singular matrix U of dot products of
   * the residual basis vectors. The solver must take a step that does
   * not use the basis vectors, and continue.
   */
   void testSingularU()
   {
      printMethod(TEST_FUNC);
      double a[4] = {1.0, 2.0, 3.0, 4.0};
      double b[4] = {2.0, 0.0, 1.0, 5.0};

      // Residuals a, a: U(0,0) = 0 for the one basis vector
      AmTestIterator iterator1;
      iterator1.residuals.resize(2);
      for (int i = 0; i < 4; ++i) {
         iterator1.residuals[0].push_back(a[i]);
         iterator1.residuals[1].push_back(a[i]);
      }
      TEST_ASSERT(iterator1.solve() == 0);
      TEST_ASSERT(iterator1.error() < 1.0E-8);

      // Residuals a, b, a: the basis vectors b - a and a - b cancel
      AmTestIterator iterator2;
      iterator2.residuals.resize(3);
      for (int i = 0; i < 4; ++i) {
         iterator2.residuals[0].push_back(a[i]);
         iterator2.residuals[1].push_back(b[i]);
         iterator2.residuals[2].push_back(a[i]);
      }
      TEST_ASSERT(iterator2.solve() == 0);
      TEST_ASSERT(iterator2.error() < 1.0E-8);

      TEST_ASSERT(amIteratorTmplTestGslErrors == 0);
   }

private:

   gsl_error_handler_t* oldHandler_;

};

TEST_BEGIN(AmIteratorTmplTest)
TEST_ADD(AmIteratorTmplTest, testSingularU)
TEST_END(AmIteratorTmplTest)

#endif
