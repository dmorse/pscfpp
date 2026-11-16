#ifndef POLYMER_CORRELATION_TEST_H
#define POLYMER_CORRELATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/tests/solvers/PolymerStub.h>

#include <pscf/correlation/Polymer.h>
#include <pscf/chem/PolymerModel.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>

#include <fstream>

using namespace Pscf;

class PolymerCorrelationTest : public UnitTest 
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

   void testReadParam(PolymerStub& p, std::string fileName) 
   {
      std::ifstream in;
      openInputFile(fileName, in);

      p.readParam(in);

      if (verbose() > 0) {
         std::cout << std::endl;
         p.writeParam(std::cout);
      }
    
      if (verbose() > 0) {
         std::cout << "\nVertices: id, size, block ids\n";
         for (int i = 0; i < p.nVertex(); ++i) {
            std::cout << i << "  " << p.vertex(i).size();
            for (int j = 0; j < p.vertex(i).size(); ++j) {
               std::cout << "  " << p.vertex(i).inPropagatorId(j)[0];
            }
            std::cout << std::endl;
         }
      }

   }

   void testDefaultConstructor()
   {
      printMethod(TEST_FUNC);
      PolymerStub p;
      Correlation::Polymer<double> c;
      c.associate(p);
   } 

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      PolymerStub p;
      Correlation::Polymer<double> c(p);
   }

   void testReadParamHomoPolymerThread() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      PolymerStub p;
      testReadParam(p, "in/HomoPolymer");
      Correlation::Polymer<double> c(p);
      int nMonomer = 2;
      c.allocate(nMonomer);
      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(length > 0.0);
      TEST_ASSERT(eq(length, 5.0));
      TEST_ASSERT(c.nBlock() == 1);
      TEST_ASSERT(c.blockIds(0).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));

      // Construct kSq array
      DArray<double> kSq;
      int nk = 10;
      kSq.allocate(nk);
      double maxQRgSq = 2.0;
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      // Allocate and zero results arrays
      DArray<double> correlations;
      DArray<double> correlationsTot;
      correlations.allocate(nk);
      correlationsTot.allocate(nk);
      for (int i = 0; i < nk; ++i) {
         correlations[i] = 0.0;
         correlationsTot[i] = 0.0;
      }

      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations);
      TEST_ASSERT(eq( correlations[0], length*length*prefactor));

      c.computeOmegaTotal(prefactor, kSq, correlationsTot);
      TEST_ASSERT(eq( correlationsTot[0], length*length*prefactor));

      double ks, x, g, value;
      for (int i = 0; i < nk; ++i) {
         ks = kSq[i];
         x = ks*RgSq;
         g = 2.0*(std::exp(-x) - 1 + x)/(x*x);
         g = prefactor*length*length*g;
         value = c.computeOmega(0, 0, prefactor, ks);
         if (i > 0) {
            TEST_ASSERT(eq(value, g));
         }
         TEST_ASSERT(eq(value, correlations[i]));
         TEST_ASSERT(eq(value, correlationsTot[i]));
         if (verbose() > 0) {
            std::cout << "\n " << Dbl(ks*RgSq) << Dbl(value);
         }
      }

   }

   void testReadParamDiblockThread() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      PolymerStub p;
      testReadParam(p, "in/PolymerDiblock");
      Correlation::Polymer<double> c(p);

      int nMonomer = 2;
      c.allocate(nMonomer);

      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(eq(length, 5.0));
      TEST_ASSERT(c.nBlock() == 2);
      TEST_ASSERT(c.blockIds(0).size() == 1);
      TEST_ASSERT(c.blockIds(1).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(c.blockIds(1)[0] == 1);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));
      TEST_ASSERT(eq(c.rSq(1,1), 0.0));
      TEST_ASSERT(eq(c.rSq(0,1), 0.0));
      TEST_ASSERT(eq(c.rSq(1,0), 0.0));

      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlationsTot;
      DArray<double> correlationsH;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlationsTot.allocate(nk);
      correlationsH.allocate(nk);
      for (int i=0; i < nk; ++i) {
         correlations00[i] = 0.0;
         correlations01[i] = 0.0;
         correlations10[i] = 0.0;
         correlations11[i] = 0.0;
         correlationsTot[i] = 0.0;
         correlationsH[i] = 0.0;
      }

      // Compute arrays of diblock copolymer properties
      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations00);
      c.computeOmega(0, 1, prefactor, kSq, correlations01);
      c.computeOmega(1, 0, prefactor, kSq, correlations10);
      c.computeOmega(1, 1, prefactor, kSq, correlations11);
      c.computeOmegaTotal(prefactor, kSq, correlationsTot);

      // Compute hompolymer properties for comparison
      PolymerStub ph;
      testReadParam(ph, "in/HomoPolymer");
      Correlation::Polymer<double> ch(ph);
      ch.allocate(nMonomer);
      ch.setup(kuhn);
      TEST_ASSERT(eq(ch.totalLength(), length));
      ch.computeOmegaTotal(prefactor, kSq, correlationsH);

      double l0, l1;
      l0 = p.block(0).length();
      l1 = p.block(1).length();
      TEST_ASSERT(eq( length, l0 + l1));
      double ks, x, c00, c01, c10, c11, cTot, sum;
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         cTot = correlationsTot[i];
         sum = c00 + c01 + c10 + c11;
         TEST_ASSERT(eq(c10, c01));
         if (i == 0) {
            TEST_ASSERT(eq(c00, prefactor*l0*l0));
            TEST_ASSERT(eq(c11, prefactor*l1*l1));
            TEST_ASSERT(eq(c01, prefactor*l0*l1));
         }
         TEST_ASSERT(eq(sum, cTot));
         TEST_ASSERT(eq(sum, correlationsH[i]));
         if (verbose() > 0) {
            ks = kSq[i];
            x = ks*RgSq;
            std::cout << "\n"  << x
                      << Dbl(c00) << Dbl(c01)
                      << Dbl(c10) << Dbl(c11)
                      << Dbl(cTot);
         }
      }

   }

   void testReadParamTriblockABAThread() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Thread);

      PolymerStub p;
      testReadParam(p, "in/PolymerTriblock");
      Correlation::Polymer<double> c(p);

      int nMonomer = 2;
      c.allocate(nMonomer);

      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(eq(length, 5.0));
      TEST_ASSERT(c.nBlock() == 3);
      TEST_ASSERT(c.blockIds(0).size() == 2);
      TEST_ASSERT(c.blockIds(1).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(c.blockIds(0)[1] == 2);
      TEST_ASSERT(c.blockIds(1)[0] == 1);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));
      TEST_ASSERT(eq(c.rSq(1,1), 0.0));
      TEST_ASSERT(eq(c.rSq(2,2), 0.0));
      TEST_ASSERT(eq(c.rSq(0,1), 0.0));
      TEST_ASSERT(eq(c.rSq(1,0), 0.0));
      TEST_ASSERT(eq(c.rSq(1,2), 0.0));
      TEST_ASSERT(eq(c.rSq(0,2), p.block(1).length()*b*b));

      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlations02;
      DArray<double> correlations12;
      DArray<double> correlations22;
      DArray<double> correlationsTot;
      DArray<double> correlationsH;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlations02.allocate(nk);
      correlations12.allocate(nk);
      correlations22.allocate(nk);
      correlationsTot.allocate(nk);
      correlationsH.allocate(nk);
      for (int i=0; i < nk; ++i) {
         correlations00[i] = 0.0;
         correlations01[i] = 0.0;
         correlations10[i] = 0.0;
         correlations11[i] = 0.0;
         correlations02[i] = 0.0;
         correlations12[i] = 0.0;
         correlations22[i] = 0.0;
         correlationsTot[i] = 0.0;
         correlationsH[i] = 0.0;
      }

      // Compute arrays of diblock copolymer properties
      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations00);
      c.computeOmega(0, 1, prefactor, kSq, correlations01);
      c.computeOmega(1, 0, prefactor, kSq, correlations10);
      c.computeOmega(1, 1, prefactor, kSq, correlations11);
      c.computeOmega(0, 2, prefactor, kSq, correlations02);
      c.computeOmega(1, 2, prefactor, kSq, correlations12);
      c.computeOmega(2, 2, prefactor, kSq, correlations22);
      c.computeOmegaTotal(prefactor, kSq, correlationsTot);

      // Compute hompolymer properties for comparison
      PolymerStub ph;
      testReadParam(ph, "in/HomoPolymer");
      Correlation::Polymer<double> ch(ph);
      ch.allocate(nMonomer);
      ch.setup(kuhn);
      TEST_ASSERT(eq(ch.totalLength(), length));
      ch.computeOmegaTotal(prefactor, kSq, correlationsH);

      double c00, c01, c10, c11, c02, c12, c22, cTot, cH, sum;
      double error1, error2, error3;
      double max1, max2, max3;
      max1 = 0.0;
      max2 = 0.0;
      max3 = 0.0;
      //setVerbose(1);
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         c02 = correlations02[i];
         c12 = correlations12[i];
         c22 = correlations22[i];
         cTot = correlationsTot[i];
         cH = correlationsH[i];
         sum = c00 + c11 + c22 + 2.0*c01 + 2.0*c12+ 2.0*c02;
         error1 = std::abs(c10 - c01);
         error2 = std::abs(cTot - cH);
         error3 = std::abs(sum - cTot);
	 if (error1 > max1) max1 = error1;
	 if (error2 > max2) max2 = error2;
	 if (error3 > max3) max3 = error3;
	 if (verbose() > 0) {
            std::cout << "\n" << Int(i) << "  " << Dbl(kSq[i]*RgSq)
		      << Dbl(sum) << Dbl(cTot) << Dbl(cH);
		      // << Dbl(error1) << Dbl(error2) << Dbl(error3);
         }
      }
      TEST_ASSERT(max1 < 1.0E-8);
      TEST_ASSERT(max2 < 1.0E-8);
      TEST_ASSERT(max3 < 1.0E-5);

   }

   void testReadParamHomoPolymerBead() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);

      PolymerStub p;
      testReadParam(p, "in/HomoPolymerBead");
      Correlation::Polymer<double> c(p);
      int nMonomer = 2;
      c.allocate(nMonomer);
      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(length > 0.0);
      TEST_ASSERT(eq(length, 10.0));
      TEST_ASSERT(c.nBlock() == 1);
      TEST_ASSERT(c.blockIds(0).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));

      // Construct kSq array
      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      // Allocate and zero all result arrays
      DArray<double> correlations;
      DArray<double> correlationsTot;
      correlations.allocate(nk);
      correlationsTot.allocate(nk);
      for (int i = 0; i < nk; ++i) {
         correlations[i] = 0.0;
         correlationsTot[i] = 0.0;
      }

      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations);
      TEST_ASSERT(eq( correlations[0], length*length*prefactor));
      c.computeOmegaTotal(prefactor, kSq, correlationsTot);

      double ks, value;
      double error1, error2;
      //setVerbose(1);
      for (int i = 0; i < nk; ++i) {
         ks = kSq[i];
         value = c.computeOmega(0, 0, prefactor, ks);
         if (verbose() > 0) {
            std::cout << "\n " << Dbl(ks*RgSq) << Dbl(correlations[i]);
         }
         error1 = std::abs(value - correlations[i]);
         error2 = std::abs(correlations[i] - correlationsTot[i]);
	 if (verbose() > 0) {
            Log::file() << "\n error1 = " << Dbl(error1);
            Log::file() << "\n error2 = " << Dbl(error2);
         }
	 TEST_ASSERT(error1 < 1.0E-6);
	 TEST_ASSERT(error2 < 1.0E-6);
      }

   }

   void testReadParamDiblockBead() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);

      PolymerStub p;
      testReadParam(p, "in/PolymerDiblockBead");
      Correlation::Polymer<double> c(p);

      int nMonomer = 2;
      c.allocate(nMonomer);

      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(eq(length, 10.0));
      TEST_ASSERT(c.nBlock() == 2);
      TEST_ASSERT(c.blockIds(0).size() == 1);
      TEST_ASSERT(c.blockIds(1).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(c.blockIds(1)[0] == 1);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));
      TEST_ASSERT(eq(c.rSq(1,1), 0.0));
      TEST_ASSERT(eq(c.rSq(0,1), b*b));
      TEST_ASSERT(eq(c.rSq(1,0), b*b));

      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlationsTot;
      DArray<double> correlationsH;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlationsTot.allocate(nk);
      correlationsH.allocate(nk);
      for (int i=0; i < nk; ++i) {
         correlations00[i] = 0.0;
         correlations01[i] = 0.0;
         correlations10[i] = 0.0;
         correlations11[i] = 0.0;
         correlationsTot[i] = 0.0;
         correlationsH[i] = 0.0;
      }

      // Compute arrays of diblock copolymer properties
      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations00);
      c.computeOmega(0, 1, prefactor, kSq, correlations01);
      c.computeOmega(1, 0, prefactor, kSq, correlations10);
      c.computeOmega(1, 1, prefactor, kSq, correlations11);
      c.computeOmegaTotal(prefactor, kSq, correlationsTot);

      // Compute hompolymer properties for comparison
      PolymerStub ph;
      testReadParam(ph, "in/HomoPolymerBead");
      Correlation::Polymer<double> ch(ph);
      ch.allocate(nMonomer);
      ch.setup(kuhn);
      TEST_ASSERT(eq(ch.totalLength(), length));
      ch.computeOmegaTotal(prefactor, kSq, correlationsH);

      double ks, x, c00, c01, c10, c11, cTot, cHT, sum;
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         cTot = correlationsTot[i];
         cHT = correlationsH[i];
         if (verbose() > 0) {
            ks = kSq[i];
            x = ks*RgSq;
            std::cout << "\n"  << x
                      << Dbl(c00) << Dbl(c01)
                      << Dbl(c10) << Dbl(c11);
         }
         TEST_ASSERT(eq(c10, c01));
         sum = c00 + c01 + c10 + c11;
         TEST_ASSERT(eq(sum, cTot));
         TEST_ASSERT(eq(sum, cHT));
      }

   }

   void testReadParamTriblockABABead() 
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);

      PolymerStub p;
      testReadParam(p, "in/PolymerTriblockBead");
      Correlation::Polymer<double> c(p);

      int nMonomer = 2;
      c.allocate(nMonomer);

      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      double b = 1.3;
      kuhn[0] = b;
      kuhn[1] = b;
      c.setup(kuhn);

      double length = c.totalLength();
      TEST_ASSERT(eq(length, 10.0));
      TEST_ASSERT(c.nBlock() == 3);
      TEST_ASSERT(c.blockIds(0).size() == 2);
      TEST_ASSERT(c.blockIds(1).size() == 1);
      TEST_ASSERT(c.blockIds(0)[0] == 0);
      TEST_ASSERT(c.blockIds(0)[1] == 2);
      TEST_ASSERT(c.blockIds(1)[0] == 1);
      TEST_ASSERT(eq(c.rSq(0,0), 0.0));
      TEST_ASSERT(eq(c.rSq(1,1), 0.0));
      TEST_ASSERT(eq(c.rSq(2,2), 0.0));
      TEST_ASSERT(eq(c.rSq(0,1), b*b));
      TEST_ASSERT(eq(c.rSq(1,0), b*b));
      TEST_ASSERT(eq(c.rSq(1,2), b*b));
      TEST_ASSERT(eq(c.rSq(0,2), (p.block(1).nBead() + 1.0)*b*b));

      // Create hompolymer object for comparison
      PolymerStub ph;
      testReadParam(ph, "in/HomoPolymerBead");
      Correlation::Polymer<double> ch(ph);
      ch.allocate(nMonomer);
      ch.setup(kuhn);
      TEST_ASSERT(eq(ch.totalLength(), length));

      // Construct kSq input array
      DArray<double> kSq;
      int nk = 10;
      double maxQRgSq = 2.0;
      kSq.allocate(nk);
      double Rg = b * sqrt(length/6.0);
      double RgSq = Rg*Rg;
      for (int i = 0; i < nk; ++i) {
         kSq[i] = maxQRgSq*double(i)/(double(nk)*RgSq);
      }

      // Allocate and zero correlation output arrays
      DArray<double> correlations00;
      DArray<double> correlations01;
      DArray<double> correlations10;
      DArray<double> correlations11;
      DArray<double> correlations02;
      DArray<double> correlations12;
      DArray<double> correlations22;
      DArray<double> correlationsTot;
      DArray<double> correlationsH;
      correlations00.allocate(nk);
      correlations01.allocate(nk);
      correlations10.allocate(nk);
      correlations11.allocate(nk);
      correlations02.allocate(nk);
      correlations12.allocate(nk);
      correlations22.allocate(nk);
      correlationsTot.allocate(nk);
      correlationsH.allocate(nk);
      for (int i=0; i < nk; ++i) {
         correlations00[i] = 0.0;
         correlations01[i] = 0.0;
         correlations10[i] = 0.0;
         correlations11[i] = 0.0;
         correlations02[i] = 0.0;
         correlations12[i] = 0.0;
         correlations22[i] = 0.0;
         correlationsTot[i] = 0.0;
         correlationsH[i] = 0.0;
      }

      // Compute arrays of diblock copolymer properties
      double prefactor = 0.7;
      c.computeOmega(0, 0, prefactor, kSq, correlations00);
      c.computeOmega(0, 1, prefactor, kSq, correlations01);
      c.computeOmega(1, 0, prefactor, kSq, correlations10);
      c.computeOmega(1, 1, prefactor, kSq, correlations11);
      c.computeOmega(0, 2, prefactor, kSq, correlations02);
      c.computeOmega(1, 2, prefactor, kSq, correlations12);
      c.computeOmega(2, 2, prefactor, kSq, correlations22);
      c.computeOmegaTotal(prefactor, kSq, correlationsTot);

      // Compute hompolymer correlation function
      ch.computeOmegaTotal(prefactor, kSq, correlationsH);

      // Compare totals
      double c00, c01, c10, c11, c02, c12, c22, cTot, cHT, sum;
      for (int i = 0; i < nk; ++i) {
         c00 = correlations00[i];
         c01 = correlations01[i];
         c10 = correlations10[i];
         c11 = correlations11[i];
         c02 = correlations02[i];
         c12 = correlations12[i];
         c22 = correlations22[i];
         cTot = correlationsTot[i];
         cHT = correlationsH[i];
         TEST_ASSERT(eq(c10, c01));
         sum = c00 + c11 + c22 + 2.0*c01 + 2.0*c12+ 2.0*c02;
         TEST_ASSERT(eq(sum, cTot));
         TEST_ASSERT(eq(sum, cHT));
      }

   }

};

TEST_BEGIN(PolymerCorrelationTest)
TEST_ADD(PolymerCorrelationTest, testDefaultConstructor)
TEST_ADD(PolymerCorrelationTest, testConstructor)
TEST_ADD(PolymerCorrelationTest, testReadParamHomoPolymerThread)
TEST_ADD(PolymerCorrelationTest, testReadParamDiblockThread)
TEST_ADD(PolymerCorrelationTest, testReadParamTriblockABAThread)
TEST_ADD(PolymerCorrelationTest, testReadParamHomoPolymerBead)
TEST_ADD(PolymerCorrelationTest, testReadParamDiblockBead)
TEST_ADD(PolymerCorrelationTest, testReadParamTriblockABABead)
TEST_END(PolymerCorrelationTest)

#endif
