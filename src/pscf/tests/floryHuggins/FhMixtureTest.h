#ifndef PSCF_FH_MIXTURE_TEST_H
#define PSCF_FH_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhMixture.h>
#include <pscf/floryHuggins/FhMolecule.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <util/containers/DArray.h>
#include <util/misc/Log.h>

#include <fstream>
#include <limits>

using namespace Pscf;
using namespace Util;

class FhMixtureTest : public UnitTest 
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
      FhMixture mixture;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixture", in);
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

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixture", in);
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

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixture", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(mixture.nMonomer());
      openInputFile("in/FhInteraction", in);
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

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixture", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(mixture.nMonomer());
      openInputFile("in/FhInteraction", in);
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

   /*
   * A/B blend of molecules of size 1 at chi = 2, its critical point.
   * At the guess phi = (0.5, 0.5), the Jacobian of the Newton steps of
   * computePhi is singular, and the GSL solver aborted the program.
   * computePhi must throw an Exception. From the guess (0.7, 0.3), it
   * must converge.
   */
   void testComputePhiSingular() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 2.0);

      DArray<double> mu, phi;
      mu.allocate(2);
      phi.allocate(2);
      mu[0] = 0.1;
      mu[1] = 0.0;
      phi[0] = 0.5;
      phi[1] = 0.5;
      double xi = 0.0;
      bool thrown = false;
      try {
         mixture.computePhi(interaction, mu, phi, xi);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      phi[0] = 0.7;
      phi[1] = 0.3;
      xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      TEST_ASSERT(eq(mixture.mu(0), mu[0]));
      TEST_ASSERT(eq(mixture.mu(1), mu[1]));
      double phiA = mixture.phi(0);
      double phiB = mixture.phi(1);
      TEST_ASSERT(phiA > 0.5);
      TEST_ASSERT(std::abs(log(phiA/phiB) + 2.0*(phiB - phiA) - 0.1)
                  < 1.0E-9);
   }

   /*
   * A/B blend of molecules of size 1 at chi = 2 with a dilute species A
   * (species 0): mu_A - mu_B = -45 gives phi_A = 3.7e-21, and -800
   * gives phi_A = exp(-802), which underflows to 0. Newton steps in phi
   * did not reach phi_A from the guess 1e-4 ("Volume fractions remain
   * out of range"), and gave an infinite Jacobian where phi underflows.
   * computePhi must converge from guesses above and below phi_A.
   */
   void testComputePhiDilute() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 2.0);

      DArray<double> mu, phi;
      mu.allocate(2);
      phi.allocate(2);
      double dMu[2] = {-45.0, -800.0};
      double guess[3] = {1.0E-4, 1.0E-30, 1.0E-300};
      double phiA, xi;
      for (int i = 0; i < 2; ++i) {
         for (int j = 0; j < 3; ++j) {
            mu[0] = dMu[i];
            mu[1] = 0.0;
            phi[0] = guess[j];
            phi[1] = 1.0 - guess[j];
            xi = 0.0;
            mixture.computePhi(interaction, mu, phi, xi);
            TEST_ASSERT(eq(mixture.mu(0), mu[0]));
            TEST_ASSERT(eq(mixture.mu(1), mu[1]));

            // phi_A = exp(mu_A - mu_B - chi), to relative order phi_A
            phiA = std::exp(dMu[i] - 2.0);
            TEST_ASSERT(std::abs(mixture.phi(0) - phiA) <= 1.0E-12*phiA);

            // f and pressure of pure B, to order phi_A ln(phi_A)
            mixture.computeFreeEnergy(interaction);
            TEST_ASSERT(eq(mixture.fHelmholtz(), -1.0));
            TEST_ASSERT(eq(mixture.pressure(), 1.0));
         }
      }
   }

   /*
   * As testComputePhiDilute, but with the dilute species last (species
   * B): the last species was the dependent one, whose phi = 1 - phi_A
   * lost the precision of phi_B, so the iteration failed. Also an ideal
   * blend (chi = 0) with mu_A - mu_B = -20 from the guess (0.6, 0.4):
   * A, the largest species of the guess, is dilute at the solution.
   */
   void testComputePhiDiluteLast() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 2.0);

      DArray<double> mu, phi;
      mu.allocate(2);
      phi.allocate(2);
      double dMu[2] = {-45.0, -800.0};
      double guess[3] = {1.0E-4, 1.0E-30, 1.0E-300};
      double phiB, xi;
      for (int i = 0; i < 2; ++i) {
         for (int j = 0; j < 3; ++j) {
            mu[0] = 0.0;
            mu[1] = dMu[i];
            phi[0] = 1.0 - guess[j];
            phi[1] = guess[j];
            xi = 0.0;
            mixture.computePhi(interaction, mu, phi, xi);
            TEST_ASSERT(eq(mixture.mu(0), mu[0]));
            TEST_ASSERT(eq(mixture.mu(1), mu[1]));
            phiB = std::exp(dMu[i] - 2.0);
            TEST_ASSERT(std::abs(mixture.phi(1) - phiB) <= 1.0E-12*phiB);
         }
      }

      // Ideal blend: phi_A/phi_B = exp(mu_A - mu_B)
      interaction.setChi(0, 1, 0.0);
      mu[0] = -20.0;
      mu[1] = 0.0;
      phi[0] = 0.6;
      phi[1] = 0.4;
      xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      TEST_ASSERT(eq(mixture.mu(0), mu[0]));
      TEST_ASSERT(eq(mixture.mu(1), mu[1]));
      double phiA = 1.0/(1.0 + std::exp(20.0));
      TEST_ASSERT(std::abs(mixture.phi(0) - phiA) <= 1.0E-12*phiA);
   }

   /*
   * Blend of A and B homopolymers and an AB diblock, all of size 1, at
   * chi = 1.5, in which the diblock is absent (a closed species with
   * phi = 0, so mu = -infinity). computeFreeEnergy gave NaN for the
   * composition (0.5, 0.5, 0), and computePhi with target mu = (0, 0,
   * -infinity) gave NaN. Both must give the A/B blend with phi = 1/2.
   */
   void testComputePhiAbsent() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixtureBlend3", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 1.5);

      DArray<double> mu, phi;
      mu.allocate(3);
      phi.allocate(3);
      double f = std::log(0.5) - 1.0 + 1.5/4.0;

      phi[0] = 0.5;
      phi[1] = 0.5;
      phi[2] = 0.0;
      mixture.setComposition(phi);
      mixture.computeMu(interaction, 0.0);
      mixture.computeFreeEnergy(interaction);
      TEST_ASSERT(eq(mixture.fHelmholtz(), f));
      TEST_ASSERT(eq(mixture.mu(0), std::log(0.5) + 0.75));
      TEST_ASSERT(eq(mixture.pressure(), -f + std::log(0.5) + 0.75));

      mu[0] = 0.0;
      mu[1] = 0.0;
      mu[2] = -std::numeric_limits<double>::infinity();
      phi[0] = 0.7;
      phi[1] = 0.3;
      double xi = 0.0;
      mixture.computePhi(interaction, mu, phi, xi);
      TEST_ASSERT(eq(mixture.phi(0), 0.5));
      TEST_ASSERT(eq(mixture.phi(1), 0.5));
      TEST_ASSERT(mixture.phi(2) == 0.0);
      TEST_ASSERT(eq(mixture.mu(0), 0.0));
      TEST_ASSERT(eq(mixture.mu(1), 0.0));
      mixture.computeFreeEnergy(interaction);
      TEST_ASSERT(eq(mixture.fHelmholtz(), f));
      TEST_ASSERT(eq(mixture.pressure(), -f));
   }

   /*
   * setComposition set the phi of the last species to 1 minus the sum
   * of the others. A dilute last species was lost (phi = 1e-20 became
   * 0, and mu = -infinity), and the check that the others add up to at
   * most 1 failed for an absent last species when the others add up
   * to 1 plus a rounding error.
   */
   void testSetCompositionDilute() {
      printMethod(TEST_FUNC);

      FhMixture mixture;
      std::ifstream in;
      openInputFile("in/FhMixtureBlend", in);
      mixture.readParam(in);
      in.close();

      FhInteraction interaction;
      interaction.setNMonomer(2);
      interaction.setChi(0, 1, 2.0);

      DArray<double> phi;
      phi.allocate(2);
      phi[0] = 1.0 - 1.0E-20;
      phi[1] = 1.0E-20;
      mixture.setComposition(phi);
      TEST_ASSERT(mixture.phi(1) == 1.0E-20);
      mixture.computeMu(interaction, 0.0);
      TEST_ASSERT(eq(mixture.mu(1), std::log(1.0E-20) + 2.0));

      FhMixture mixture3;
      openInputFile("in/FhMixtureBlend3", in);
      mixture3.readParam(in);
      in.close();
      interaction.setChi(0, 1, 1.5);

      // phi3[0] + phi3[1] = 1 + 2.2e-16
      DArray<double> phi3;
      phi3.allocate(3);
      phi3[0] = 0.5;
      phi3[1] = std::nextafter(0.5, 1.0);
      phi3[2] = 0.0;
      mixture3.setComposition(phi3);
      TEST_ASSERT(mixture3.phi(2) == 0.0);
      mixture3.computeMu(interaction, 0.0);
      mixture3.computeFreeEnergy(interaction);
      TEST_ASSERT(eq(mixture3.fHelmholtz(), std::log(0.5) - 1.0 + 1.5/4.0));
   }

};

TEST_BEGIN(FhMixtureTest)
TEST_ADD(FhMixtureTest, testConstructor)
TEST_ADD(FhMixtureTest, testReadWrite)
TEST_ADD(FhMixtureTest, testSetComposition)
TEST_ADD(FhMixtureTest, testComputeMu)
TEST_ADD(FhMixtureTest, testComputePhi)
TEST_ADD(FhMixtureTest, testComputePhiSingular)
TEST_ADD(FhMixtureTest, testComputePhiDilute)
TEST_ADD(FhMixtureTest, testComputePhiDiluteLast)
TEST_ADD(FhMixtureTest, testComputePhiAbsent)
TEST_ADD(FhMixtureTest, testSetCompositionDilute)
TEST_END(FhMixtureTest)

#endif
