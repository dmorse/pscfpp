#ifndef PSPC_FFTW_TEST_H
#define PSPC_FFTW_TEST_H

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
//using namespace Pspc;

class cuFFTTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
//    void testFFTW_1D() {
//       printMethod(TEST_FUNC);
//       printEndl();

//       DArray<double> in;
//       DArray<fftw_complex> out;
//       int n = 10;
//       in.allocate(n);
//       out.allocate(n/2+1);
//       unsigned int flags = FFTW_ESTIMATE;

//       double x;
//       double twoPi = 2.0*Constants::Pi;
//       for (int i = 0; i < n; ++i) {
//          x = twoPi*float(i)/float(n); 
//          in[i] = cos(x);
//          std::cout << Dbl(in[i]);
//       }
//       std::cout << std::endl;

//       //std::cout << "Entering forward transform plan creation" << std::endl;
//       fftw_plan plan_f = fftw_plan_dft_r2c_1d(n, &in[0], &out[0], flags);
//       //std::cout << "Finished forward plan creation" << std::endl;
//       fftw_execute(plan_f);
//       //std::cout << "Finished forward transform" << std::endl;
       
//       fftw_plan plan_r = fftw_plan_dft_c2r_1d(n, &out[0], &in[0], flags);
//       //std::cout << "Finished inverse plan creation" << std::endl;
//       fftw_execute(plan_r);
//       //std::cout << "Finished inverse transform" << std::endl;
    
//       double factor = 1.0/double(n); 
//       for (int i = 0; i < n/2 + 1; ++i) {
//          out[i][0] *= factor;
//          out[i][1] *= factor;
//          //std::cout << out[i][0] << "  " << out[i][1] << std::endl;
//       }

//       for (int i = 0; i < n; ++i) {
//          std::cout << Dbl(in[i]);
//       }
//       std::cout << std::endl;
//    }

// };

TEST_BEGIN(cuFFTTest)
// TEST_ADD(cuFFTTest, testFFTW_1D)
TEST_END(cuFFTTest)

#endif
