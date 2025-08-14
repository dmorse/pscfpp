#ifndef PSCF_FLORY_HUGGINS_MIXTURE_TEST_H
#define PSCF_FLORY_HUGGINS_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/Mixture.h>
#include <pscf/floryHuggins/Molecule.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/DArray.h>
#include <util/misc/Log.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class MixtureTest : public UnitTest 
{

public:

   void setUp()
   {
      //setVerbose(1);
   }

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      FloryHuggins::Mixture mixture;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);

      FloryHuggins::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);

      TEST_ASSERT(mixture.nMolecule() == 2);

      TEST_ASSERT(eq(mixture.molecule(0).size(), 5.0));
      TEST_ASSERT(mixture.molecule(0).nClump() == 2);
      TEST_ASSERT(mixture.molecule(0).clump(0).monomerId() == 0);
      TEST_ASSERT(eq(mixture.molecule(0).clump(0).size(), 2.0));
      TEST_ASSERT(mixture.molecule(0).clump(1).monomerId() == 1);
      TEST_ASSERT(eq(mixture.molecule(0).clump(1).size(), 3.0));

      TEST_ASSERT(eq(mixture.molecule(1).size(), 1.0));
      TEST_ASSERT(mixture.molecule(1).nClump() == 1);
      TEST_ASSERT(mixture.molecule(1).clump(0).monomerId() == 0);
      TEST_ASSERT(eq(mixture.molecule(1).clump(0).size(), 1.0));

      if (verbose() > 0) {
         printEndl();
         mixture.writeParam(Log::file());
      }
   }

   void testSetComposition() {
      printMethod(TEST_FUNC);

      FloryHuggins::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);

      DArray<double> phi;
      phi.allocate(2);
      phi[0] = 0.6;
      phi[1] = 0.4;
      mixture.setComposition(phi);

      TEST_ASSERT(eq(mixture.phi(0), 0.6));
      TEST_ASSERT(eq(mixture.phi(1), 0.4));
      TEST_ASSERT(eq(mixture.c(0), 0.64));
      TEST_ASSERT(eq(mixture.c(1), 0.36));

   }

   void testComputeMu() {
      printMethod(TEST_FUNC);

      FloryHuggins::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);
      in.close();

      Interaction interaction;
      interaction.setNMonomer(mixture.nMonomer());
      openInputFile("in/Interaction", in);
      interaction.readParam(in);
      in.close();

      DArray<double> phi;
      phi.allocate(2);
      phi[0] = 0.6;
      phi[1] = 0.4;
      mixture.setComposition(phi);

      double xi = 3.0;
      mixture.computeMu(interaction, xi);

      TEST_ASSERT(eq(mixture.phi(0), 0.6));
      TEST_ASSERT(eq(mixture.phi(1), 0.4));
      TEST_ASSERT(eq(mixture.c(0), 0.64));
      TEST_ASSERT(eq(mixture.c(1), 0.36));

      double mu0, mu1, w0, w1, chi;
      chi = 2.0;
      w0 = mixture.c(1)*chi;
      w1 = mixture.c(0)*chi;
      mu0 = log(0.6) + w0*2.0 + w1*3.0 + xi*5.0;
      TEST_ASSERT(eq(mixture.mu(0), mu0));
      mu1 = log(0.4) + w0*1.0 + xi*1.0;
      TEST_ASSERT(eq(mixture.mu(1), mu1));

   }

   void testComputePhi() {
      printMethod(TEST_FUNC);

      FloryHuggins::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);
      in.close();

      Interaction interaction;
      interaction.setNMonomer(mixture.nMonomer());
      openInputFile("in/Interaction", in);
      interaction.readParam(in);
      in.close();

      DArray<double> phi;
      phi.allocate(2);
      phi[0] = 0.6;
      phi[1] = 0.4;
      mixture.setComposition(phi);

      double xi = 3.0;
      mixture.computeMu(interaction, xi);

      TEST_ASSERT(eq(mixture.phi(0), 0.6));
      TEST_ASSERT(eq(mixture.phi(1), 0.4));
      TEST_ASSERT(eq(mixture.c(0), 0.64));
      TEST_ASSERT(eq(mixture.c(1), 0.36));

      DArray<double> mu;
      mu.allocate(2);
      mu[0] = mixture.mu(0) + 0.10;
      mu[1] = mixture.mu(1) - 0.05;
      xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);

      // Note: Throw exception if convergence fails, so
      // normal completion already indicates success.
 
      TEST_ASSERT(eq(mixture.mu(0), mu[0]));
      TEST_ASSERT(eq(mixture.mu(1), mu[1]));

      mixture.computeFreeEnergy(interaction);
      if (verbose() > 0) {
         printEndl();
         Log::file() << "fHelmholtz = " << mixture.fHelmholtz() << "\n";
         Log::file() << "pressure   = " << mixture.pressure() << "\n";
      }
   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor)
TEST_ADD(MixtureTest, testReadWrite)
TEST_ADD(MixtureTest, testSetComposition)
TEST_ADD(MixtureTest, testComputeMu)
TEST_ADD(MixtureTest, testComputePhi)
TEST_END(MixtureTest)

#endif
