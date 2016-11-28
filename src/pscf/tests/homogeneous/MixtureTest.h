#ifndef PSCF_HOMOGENEOUS_MIXTURE_TEST_H
#define PSCF_HOMOGENEOUS_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/homogeneous/Mixture.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/containers/DArray.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

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
      Homogeneous::Mixture mixture;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Homogeneous::Mixture mixture;
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
      mixture.writeParam(std::cout) ;
   }

   void testSetComposition() {
      printMethod(TEST_FUNC);
      printEndl();

      Homogeneous::Mixture mixture;
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
      printEndl();

      Homogeneous::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);
      in.close();

      ChiInteraction interaction;
      interaction.setNMonomer(mixture.nMonomer());
      openInputFile("in/ChiInteraction", in);
      interaction.readParam(in);
      in.close();

      DArray<double> phi;
      phi.allocate(2);
      phi[0] = 0.6;
      phi[1] = 0.4;
      double xi = 3.0;
      mixture.computeMu(interaction, phi, xi);

      TEST_ASSERT(eq(mixture.phi(0), 0.6));
      TEST_ASSERT(eq(mixture.phi(1), 0.4));
      TEST_ASSERT(eq(mixture.c(0), 0.64));
      TEST_ASSERT(eq(mixture.c(1), 0.36));

      double mu0, mu1, w0, w1, chi;
      chi = 2.0;
      w0 = mixture.c(1)*chi;
      w1 = mixture.c(0)*chi;
      mu0 = log(0.6) + w0*2.0 + w1*3.0 + xi*5.0;
      //std::cout << "Mu(0) = " << mixture.mu(0) << std::endl;
      //std::cout << "Mu(0) = " << mu0 << std::endl;
      TEST_ASSERT(eq(mixture.mu(0), mu0));
      mu1 = log(0.4) + w0*1.0 + xi*1.0;
      //std::cout << "Mu(1) = " << mixture.mu(1) << std::endl;
      //std::cout << "Mu(1) = " << mu1 << std::endl;
      TEST_ASSERT(eq(mixture.mu(1), mu1));

   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor)
TEST_ADD(MixtureTest, testReadWrite)
TEST_ADD(MixtureTest, testSetComposition)
TEST_ADD(MixtureTest, testComputeMu)
TEST_END(MixtureTest)

#endif
