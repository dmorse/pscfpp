#ifndef MIXTURE_CORRELATION_TEST_H
#define MIXTURE_CORRELATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/tests/solvers/MixtureStub.h>
#include <pscf/tests/solvers/PolymerStub.h>

#include <pscf/correlation/Mixture.h>
#include <pscf/correlation/Polymer.h>
#include <pscf/chem/PolymerModel.h>
#include <util/format/Dbl.h>

#include <fstream>

using namespace Pscf;

class MixtureCorrelationTest : public UnitTest 
{

public:

   void setUp()
   {  
      setVerbose(0); 
   }

   void tearDown()
   {  
      setVerbose(0); 
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testReadParam(MixtureStub& p, std::string fileName) 
   {
      std::ifstream in;
      openInputFile(fileName, in);

      p.readParam(in);

      if (verbose() > 0) {
         std::cout << std::endl;
         p.writeParam(std::cout);
      }

   }

   void testDefaultConstructor()
   {
      printMethod(TEST_FUNC);
      MixtureStub p;
      Correlation::Mixture<double> c;
      c.associate(p);
      TEST_ASSERT(!c.isAllocated());
   } 

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      MixtureStub p;
      Correlation::Mixture<double> c(p);
      TEST_ASSERT(!c.isAllocated());
   }


   void testReadParamDiblockThread() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      MixtureStub mixture;
      testReadParam(mixture, "in/MixtureDiblockThread");
      Correlation::Mixture<double> c(mixture);
      c.allocate();
      c.setup();
      TEST_ASSERT(c.isAllocated());

      // Test: Mixture of a homopolymer and diblock of equal length 5.0
      double length = c.polymer(0).totalLength();
      TEST_ASSERT(eq(length, 5.0));
      length = c.polymer(1).totalLength();
      TEST_ASSERT(eq(length, 5.0));

      // Homopolymer
      TEST_ASSERT(c.polymer(0).nBlock() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0)[0] == 0);

      // Diblock
      TEST_ASSERT(c.polymer(1).nBlock() == 2);
      TEST_ASSERT(c.polymer(1).blockIds(0).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(1).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0)[0] == 0);
      TEST_ASSERT(c.polymer(1).blockIds(1)[0] == 1);

      // Test for equal statistical segment lengths
      Composition<double> const & mixc = mixture;
      double b = mixc.monomer(0).kuhn();
      TEST_ASSERT(eq(mixc.monomer(1).kuhn(), b));

      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      // Allocate correlation arrays for monomer type pairs and total
      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlationsTot;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlationsTot.allocate(nk);

      // Compute arrays of diblock copolymer properties
      c.computeOmega(0, 0, kSq, correlations00);
      c.computeOmega(0, 1, kSq, correlations01);
      c.computeOmega(1, 0, kSq, correlations10);
      c.computeOmega(1, 1, kSq, correlations11);
      c.computeOmegaTotal(kSq, correlationsTot);

      // Allocate, setup and compute properties for homopolymer melt
      MixtureStub mixtureH;
      testReadParam(mixtureH, "in/MixtureHomopolymerThread");
      Correlation::Mixture<double> cH(mixture);
      cH.allocate();
      cH.setup();
      TEST_ASSERT(eq(cH.polymer(0).totalLength(), length));
      DArray<double> correlationsH;
      correlationsH.allocate(nk);
      cH.computeOmegaTotal(kSq, correlationsH);

      double vMonomer = mixture.vMonomer();
      double phi0 = mixture.polymer(0).phi();
      double phi1 = mixture.polymer(1).phi();
      TEST_ASSERT(eq(1.0, phi0 + phi1));

      // Loop over array results
      double ks, x, c00, c01, c10, c11, cTot, sum, cHT, g;
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         cTot = correlationsTot[i];
         cHT = correlationsH[i];
         TEST_ASSERT(eq(c10, c01));
         sum = c00 + c01 + c10 + c11;
         TEST_ASSERT(eq(sum, cTot));
         TEST_ASSERT(eq(sum, cHT));
         ks = kSq[i];
         x = ks*RgSq;
         if (i == 0) {
            TEST_ASSERT(eq(sum, length/vMonomer));
         } else {
            g = 2.0*(std::exp(-x) - 1.0 + x)/(x*x);
            g = length*g/vMonomer;
            TEST_ASSERT(eq(sum, g));
         }
         if (verbose() > 0) {
            std::cout << "\n"  << x
                      << Dbl(c00) << Dbl(c01)
                      << Dbl(c10) << Dbl(c11)
                      << Dbl(cTot) << Dbl(cHT);
         }
      }

   }

   void testReadParamDiblockBead() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);

      MixtureStub mixture;
      testReadParam(mixture, "in/MixtureDiblockBead");
      Correlation::Mixture<double> c(mixture);
      TEST_ASSERT(!c.isAllocated());
      c.allocate();
      c.setup();

      // Test: Mixture of a homopolymer and diblock of equal length 10
      TEST_ASSERT(c.isAllocated());
      double length = c.polymer(0).totalLength();
      TEST_ASSERT(eq(length, 10.0));
      length = c.polymer(1).totalLength();
      TEST_ASSERT(eq(length, 10.0));

      // Homopolymer
      TEST_ASSERT(c.polymer(0).nBlock() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0)[0] == 0);

      // Diblock
      TEST_ASSERT(c.polymer(1).nBlock() == 2);
      TEST_ASSERT(c.polymer(1).blockIds(0).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(1).size() == 1);
      TEST_ASSERT(c.polymer(1).blockIds(0)[0] == 0);
      TEST_ASSERT(c.polymer(1).blockIds(1)[0] == 1);

      // Test for equal statistical segment lengths
      Composition<double> const & mixc = mixture;
      double b = mixc.monomer(0).kuhn();
      TEST_ASSERT(eq(mixc.monomer(1).kuhn(), b));

      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      // Allocate correlation arrays for monomer type pairs and total
      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlationsTot;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlationsTot.allocate(nk);

      // Compute arrays of diblock copolymer properties
      c.computeOmega(0, 0, kSq, correlations00);
      c.computeOmega(0, 1, kSq, correlations01);
      c.computeOmega(1, 0, kSq, correlations10);
      c.computeOmega(1, 1, kSq, correlations11);
      c.computeOmegaTotal(kSq, correlationsTot);

      // Allocate, setup and compute properties for homopolymer melt
      MixtureStub mixtureH;
      testReadParam(mixtureH, "in/MixtureHomopolymerBead");
      Correlation::Mixture<double> cH(mixture);
      cH.allocate();
      cH.setup();
      TEST_ASSERT(eq(cH.polymer(0).totalLength(), length));
      DArray<double> correlationsH;
      correlationsH.allocate(nk);
      cH.computeOmegaTotal(kSq, correlationsH);

      // Check volume fractions
      double vMonomer = mixture.vMonomer();
      double phi0 = mixture.polymer(0).phi();
      double phi1 = mixture.polymer(1).phi();
      TEST_ASSERT(eq(1.0, phi0 + phi1));
      TEST_ASSERT(eq(mixtureH.polymer(0).phi(), 1.0));

      // Loop over array results
      double ks, x, c00, c01, c10, c11, cTot, sum, cHT;
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         cTot = correlationsTot[i];
         cHT = correlationsH[i];
         TEST_ASSERT(eq(c10, c01));
         sum = c00 + c01 + c10 + c11;
         TEST_ASSERT(eq(sum, cTot));
         TEST_ASSERT(eq(sum, cHT));
         ks = kSq[i];
         x = ks*RgSq;
         if (i == 0) {
            TEST_ASSERT(eq(sum, length/vMonomer));
         }
         if (verbose() > 0) {
            std::cout << "\n"  << x
                      << Dbl(c00) << Dbl(c01)
                      << Dbl(c10) << Dbl(c11)
                      << Dbl(cTot) << Dbl(cHT);
         }
      }

   }
};

TEST_BEGIN(MixtureCorrelationTest)
TEST_ADD(MixtureCorrelationTest, testDefaultConstructor)
TEST_ADD(MixtureCorrelationTest, testConstructor)
TEST_ADD(MixtureCorrelationTest, testReadParamDiblockThread)
TEST_ADD(MixtureCorrelationTest, testReadParamDiblockBead)
TEST_END(MixtureCorrelationTest)

#endif
