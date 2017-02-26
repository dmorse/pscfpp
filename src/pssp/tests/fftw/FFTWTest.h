#ifndef PSSP_FFTW_TEST_H
#define PSSP_FFTW_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

#include <fftw3.h>

#include <iostream>
#include <fstream>

using namespace Util;
//using namespace Pscf;
//using namespace Pssp;

class FFTWTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   void testFFTW_1D() {
      printMethod(TEST_FUNC);
      printEndl();

      DArray<double> in;
      DArray<fftw_complex> out;
      int n = 10;
      in.allocate(n);
      out.allocate(n/2+1);
      unsigned int flags = FFTW_ESTIMATE;

      double x;
      double twoPi = 2.0*Constants::Pi;
      for (int i = 0; i < n; ++i) {
         x = twoPi*float(i)/float(n); 
         in[i] = cos(x);
         std::cout << Dbl(in[i]);
      }
      std::cout << std::endl;

      fftw_plan plan_f = fftw_plan_dft_r2c_1d(n, &in[0], &out[0], flags);
      fftw_execute(plan_f);
      //fftw_plan plan = fftw_plan_dft_r2c_1d(n, &in[0], &out[0], flags);
       
      fftw_plan plan_r = fftw_plan_dft_c2r_1d(n, &out[0], &in[0], flags);
      fftw_execute(plan_r);
     
      double factor = 1.0/double(n); 
      for (int i = 0; i < n; ++i) {
         out[i][0] *= factor;
         out[i][1] *= factor;
      }

      for (int i = 0; i < n; ++i) {
         std::cout << Dbl(in[i]);
      }
      std::cout << std::endl;
   }

};

TEST_BEGIN(FFTWTest)
TEST_ADD(FFTWTest, testFFTW_1D)
TEST_END(FFTWTest)

#endif
