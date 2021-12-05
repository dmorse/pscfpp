#ifndef PSPC_SWEEP_TEST_H
#define PSPC_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/SweepFactory.h>
#include <pspc/sweep/LinearSweep.h>
#include <pspc/field/BFieldComparison.h>
#include <util/format/Dbl.h>

#include <fstream>
#include <sstream>


using namespace Util;
using namespace Pscf;
using namespace Pspc;

class SweepTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   void testConstructors()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      LinearSweep<3> ls(system);
      SweepFactory<3> sf(system);
   }

   void testFactory() 
   {
      printMethod(TEST_FUNC);

      System<3> system;
      Sweep<3>* sweepPtr;
      SweepFactory<3> sf(system);
      
      sweepPtr = sf.factory("LinearSweep");
   }

   void testParameterRead() 
   {
      printMethod(TEST_FUNC);

      // Set up system with some data
      System<1> system;
      SweepTest::SetUpSystem(system, "in/block/param");
      // Set up SweepParameter objects 
      DArray< SweepParameter<1> > ps;
      ps.allocate(4);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
      }
      // Open test input file
      std::ifstream in;

      // Read in data
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         in >> ps[i];
      }

      // Assert that it is read correctly
      TEST_ASSERT(ps[0].type()=="block_length");
      TEST_ASSERT(ps[0].id(0)==0);
      TEST_ASSERT(ps[0].id(1)==0);
      TEST_ASSERT(ps[0].change()==0.25);
      TEST_ASSERT(ps[1].type()=="chi");
      TEST_ASSERT(ps[1].id(0)==0);
      TEST_ASSERT(ps[1].id(1)==1);
      TEST_ASSERT(ps[1].change()==5.00);
      TEST_ASSERT(ps[2].type()=="kuhn");
      TEST_ASSERT(ps[2].id(0)==0);
      TEST_ASSERT(ps[2].change()==0.1);
      TEST_ASSERT(ps[3].type()=="phi_polymer");
      TEST_ASSERT(ps[3].id(0)==0);
      TEST_ASSERT(ps[3].change()==-0.01);
   }

   void testParameterGet()
   {
      printMethod(TEST_FUNC);

      // Set up system
      System<1> system;
      SweepTest::SetUpSystem(system, "in/block/param");

      // Set up SweepParameter objects 
      DArray< SweepParameter<1> > ps;
      ps.allocate(4);
      std::ifstream in;
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
         in >> ps[i];
      }

      DArray<double> sysval, paramval;
      sysval.allocate(4);
      paramval.allocate(4);

      // Call get_ function to get value through parameter
      for (int i = 0; i < 4; ++i) {
         paramval[i] = ps[i].current();
      }

      // Manually check equality for each one
      sysval[0] = system.mixture().polymer(0).block(0).length();
      sysval[1] = system.interaction().chi(0,1);
      sysval[2] = system.mixture().monomer(0).step();
      sysval[3] = system.mixture().polymer(0).phi();
      for (int i = 0; i < 4; ++i) {
         TEST_ASSERT(sysval[i] == paramval[i]);
      }
      
   }

   void testParameterSet()
   {
      printMethod(TEST_FUNC);

      // Set up system
      System<1> system;
      SweepTest::SetUpSystem(system, "in/block/param");

      // Set up SweepParameter objects 
      DArray< SweepParameter<1> > ps;
      ps.allocate(4);
      std::ifstream in;
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
         in >> ps[i];
         ps[i].getInitial();
      }

      DArray<double> sysval, paramval;
      sysval.allocate(4);
      paramval.allocate(4);

      // Set for some arbitrary value of s in [0,1]
      double s = 0.295586;
      double newVal;
      for (int i = 0; i < 4; ++i) {
         newVal = ps[i].initial() + s*ps[i].change();
         ps[i].update(newVal);
      }
      // Calculate expected value of parameter using s
      for (int i = 0; i < 4; ++i) {
         paramval[i] = ps[i].initial() + s*ps[i].change();
      }

      // Manually check equality for each one
      sysval[0] = system.mixture().polymer(0).block(0).length();
      sysval[1] = system.interaction().chi(0,1);
      sysval[2] = system.mixture().monomer(0).step();
      sysval[3] = system.mixture().polymer(0).phi();
      for (int i = 0; i < 4; ++i) {
         TEST_ASSERT(sysval[i]==paramval[i]);
      }
   }

   void testLinearSweepRead()
   {
      printMethod(TEST_FUNC);

      // Set up system with Linear Sweep Object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/block/param");
   }

   void testLinearSweepBlock()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepBlock");

      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/block/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/block/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/block/w.bf");
      system.sweep();

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/block/" + std::to_string(i) +"_w.bf");
      }

      // Compare output
      DArray< BFieldComparison > comparisons;
      comparisons.allocate(5);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparisons[i].compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparisons[i].maxDiff() > maxDiff) {
            maxDiff = comparisons[i].maxDiff();
         }
      }
      setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "maxDiff = " << Dbl(maxDiff, 14, 6) << std::endl;
      }
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepChi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepChi");

      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/chi/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/chi/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/chi/w.bf");
      system.sweep();

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/chi/" + std::to_string(i) +"_w.bf");
      }

      // Compare output
      DArray< BFieldComparison > comparisons;
      comparisons.allocate(5);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparisons[i].compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparisons[i].maxDiff() > maxDiff) {
            maxDiff = comparisons[i].maxDiff();
         }
      }
      setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "maxDiff = " << Dbl(maxDiff,14,6) << std::endl;
      }
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepKuhn()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepKuhn");

      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/kuhn/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/kuhn/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/kuhn/w.bf");
      system.sweep();

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/kuhn/" + std::to_string(i) +"_w.bf");
      }

      // Compare output
      DArray< BFieldComparison > comparisons;
      comparisons.allocate(5);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparisons[i].compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparisons[i].maxDiff() > maxDiff) {
            maxDiff = comparisons[i].maxDiff();
         }
      }
      setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "maxDiff = " << Dbl(maxDiff,14,6) << std::endl;
      }
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepPhi()   
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepPhi");

      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/phi/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/phi/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/phi/w.bf");
      system.sweep();

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/phi/" + std::to_string(i) +"_w.bf");
      }

      // compare output
      DArray< BFieldComparison > comparisons;
      comparisons.allocate(5);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparisons[i].compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparisons[i].maxDiff() > maxDiff) {
            maxDiff = comparisons[i].maxDiff();
         }
      }
      setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "maxDiff = " << Dbl(maxDiff, 14, 6) << std::endl;
      }
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepSolvent()   
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepSolvent");

      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/solvent/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/solvent/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/solvent/w.bf");
      system.sweep();

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/solvent/" + std::to_string(i) +"_w.bf");
      }

      // Compare output
      DArray< BFieldComparison > comparisons;
      comparisons.allocate(5);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparisons[i].compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparisons[i].maxDiff() > maxDiff) {
            maxDiff = comparisons[i].maxDiff();
         }
      }
      setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "maxDiff = " << Dbl(maxDiff, 14, 6) << std::endl;
      }
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void SetUpSystem(System<1>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructors)
TEST_ADD(SweepTest, testFactory)
TEST_ADD(SweepTest, testParameterRead)
TEST_ADD(SweepTest, testParameterGet)
TEST_ADD(SweepTest, testParameterSet)
TEST_ADD(SweepTest, testLinearSweepRead)
TEST_ADD(SweepTest, testLinearSweepBlock)
TEST_ADD(SweepTest, testLinearSweepChi)
TEST_ADD(SweepTest, testLinearSweepKuhn)
TEST_ADD(SweepTest, testLinearSweepPhi)
TEST_ADD(SweepTest, testLinearSweepSolvent)
TEST_END(SweepTest)

#endif
