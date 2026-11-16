#ifndef R1_MIXTURE_TEST_H
#define R1_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/solvers/Solvent.h>
#include <r1/solvers/Propagator.h>
#include <util/tests/LogFileUnitTest.h>

#include <cmath>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

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

   /*
   * Read a mixture and its domain, and allocate arrays of w and c fields.
   */
   void setupMixture(char const * filename, Mixture& mix, Domain& domain,
                     DArray<Mixture::FieldT>& wFields,
                     DArray<Mixture::FieldT>& cFields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      mix.readParam(in);
      domain.readParam(in);
      in.close();
      mix.setDomain(domain);
      int nMonomer = mix.nMonomer();
      int nx = domain.nx();
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(nx);
         cFields[i].allocate(nx);
      }
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

   /*
   * Solvents of size 1800 in a uniform field 0.45, whose partition
   * function q = exp(-810) is below the range of a double: an open
   * solvent of monomer type 0 with mu = 809.5, for which exp(mu) is
   * above that range, and a closed solvent of monomer type 1 with
   * phi = 0.4. Then the closed solvent in a field alternating between
   * -0.05 and 0.95, which changes its Boltzmann weight by exp(1800).
   */
   void testSolveLargeSolvent()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/MixtureTestSolveLargeSolvent.log");

      Mixture mix;
      Domain domain;
      DArray<Mixture::FieldT> wFields;
      DArray<Mixture::FieldT> cFields;
      setupMixture("in/MixtureSolvent", mix, domain, wFields, cFields);
      int nx = domain.nx();

      for (int i = 0; i < nx; ++i) {
         wFields[0][i] = 0.45;
         wFields[1][i] = 0.45;
      }
      mix.compute(wFields, cFields);
      Solvent const & open = mix.solvent(0);
      Solvent const & closed = mix.solvent(1);
      TEST_ASSERT(std::abs(open.phi() - std::exp(-0.5)) < 1.0E-10);
      TEST_ASSERT(std::abs(closed.mu() - std::log(0.4) - 810.0) < 1.0E-10);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::abs(cFields[0][i] - 0.1 - std::exp(-0.5))
                     < 1.0E-10);
         TEST_ASSERT(std::abs(cFields[1][i] - 0.4) < 1.0E-10);
      }

      // The points with the field -0.05 hold half of the volume (with
      // half weights at the two end points of the grid): q = exp(90)/2
      for (int i = 0; i < nx; ++i) {
         wFields[1][i] = (i % 2 == 0) ? -0.05 : 0.95;
      }
      mix.compute(wFields, cFields);
      TEST_ASSERT(std::abs(closed.mu() - std::log(0.4) + 90.0
                           - std::log(2.0)) < 1.0E-10);
      for (int i = 0; i < nx; ++i) {
         double c = (i % 2 == 0) ? 0.8 : 0.0;
         TEST_ASSERT(std::abs(cFields[1][i] - c) < 1.0E-10);
      }
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
      DArray<Mixture::FieldT> wFields;
      DArray<Mixture::FieldT> cFields;
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

      double q00, q01, q10, q11;
      mix.polymer(0).propagator(0, 0).computeQ(q00);
      mix.polymer(0).propagator(0, 1).computeQ(q01);
      mix.polymer(0).propagator(1, 0).computeQ(q10);
      mix.polymer(0).propagator(1, 1).computeQ(q11);
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
TEST_ADD(MixtureTest, testSolveLargeSolvent)
TEST_ADD(MixtureTest, testSolve)
TEST_END(MixtureTest)

#endif
