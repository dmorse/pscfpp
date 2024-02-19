#ifndef R1D_MIXTURE_TEST_H
#define R1D_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/solvers/Mixture.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1d;

class MixtureTest : public LogFileUnitTest 
{

public:

   void setUp()
   {  ParamComponent::setEcho(true); }

   void tearDown()
   {
      closeLogFile();
      ParamComponent::setEcho(false);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Mixture mix;
   }

   void testReadParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MixtureTestReadParameters.log");

      std::ifstream in;
      openInputFile("in/Mixture", in);

      Mixture mix;
      mix.readParam(in);

      Domain domain;
      domain.readParam(in);
      mix.setDomain(domain);

      Log::file() << "\n";
      mix.writeParam(Log::file());
      domain.writeParam(Log::file());
   }

   void testReadParameters2()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MixtureTestReadParameters2.log");

      std::ifstream in;
      openInputFile("in/Mixture2", in);

      Mixture mix;
      mix.readParam(in);

      Domain domain;
      domain.readParam(in);
      mix.setDomain(domain);

      TEST_ASSERT(eq(mix.vMonomer(), 0.05));

      Log::file() << "\n";
      mix.writeParam(Log::file());
      domain.writeParam(Log::file());
   }

   void testSolve()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MixtureTestSolve.log");

      std::ifstream in;
      openInputFile("in/Mixture", in);

      Mixture mix;
      Domain domain;
      mix.readParam(in);
      domain.readParam(in);
      mix.setDomain(domain);

      Log::file() << "\n";
      mix.writeParam(Log::file());
      domain.writeParam(Log::file());

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

      double q00 = mix.polymer(0).propagator(0, 0).computeQ();
      double q01 = mix.polymer(0).propagator(0, 1).computeQ();
      double q10 = mix.polymer(0).propagator(1, 0).computeQ();
      double q11 = mix.polymer(0).propagator(1, 1).computeQ();
      TEST_ASSERT(abs(q01 - q00) < 1.0E-5);
      TEST_ASSERT(abs(q10 - q00) < 1.0E-5);
      TEST_ASSERT(abs(q11 - q00) < 1.0E-5);

      // Test if same Q is obtained from different methods
      Log::file() << "Propagator(0,0), Q = " << q00 << "\n";
      Log::file() << "Propagator(0,1), Q = " << q10 << "\n";
      Log::file() << "Propagator(1,0), Q = " << q01 << "\n";
      Log::file() << "Propagator(0,1), Q = " << q11 << "\n";

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
      
   }
};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor)
TEST_ADD(MixtureTest, testReadParameters)
TEST_ADD(MixtureTest, testReadParameters2)
TEST_ADD(MixtureTest, testSolve)
TEST_END(MixtureTest)

#endif
