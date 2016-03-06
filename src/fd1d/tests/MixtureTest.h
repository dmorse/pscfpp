#ifndef FD1D_MIXTURE_TEST_H
#define FD1D_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/Mixture.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Fd1d;

class MixtureTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Mixture mix;
   }

   void testReadParameters()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/Mixture", in);

      Mixture mix;
      mix.readParam(in);

      Domain domain;
      domain.readParam(in);
      mix.setDomain(domain);

      std::cout << "\n";
      mix.writeParam(std::cout);
      domain.writeParam(std::cout);
   }

   void testSolve()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/Mixture", in);

      Mixture mix;
      Domain domain;
      mix.readParam(in);
      domain.readParam(in);
      mix.setDomain(domain);

      std::cout << "\n";
      mix.writeParam(std::cout);
      domain.writeParam(std::cout);

      int nMonomer = mix.nMonomer();
      DArray<Mixture::WField> wFields;
      DArray<Mixture::CField> cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)domain.nx();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(nx);
         cFields[i].allocate(nx);
      }

      double cs;
      for (int i = 0; i < nx; ++i) {
         //cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }
      mix.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      std::cout << "Propagator(0,0), Q = " 
                << mix.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << mix.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << mix.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << mix.polymer(0).propagator(0, 1).computeQ() << "\n";

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      
   }
};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor)
TEST_ADD(MixtureTest, testReadParameters)
TEST_ADD(MixtureTest, testSolve)
TEST_END(MixtureTest)

#endif
