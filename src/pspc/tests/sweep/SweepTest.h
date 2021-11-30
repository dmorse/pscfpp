#ifndef PSPC_SWEEP_TEST_H
#define PSPC_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/SweepFactory.h>
#include <pspc/sweep/LinearSweep.h>

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
   {}

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

      // set up system with some data
      System<3> system;
      SweepTest::SetUpSystem(system, "in/diblock/bcc/param.flex");
      // set up LinearSweepParameter objects 
      DArray< LinearSweepParameter<3> > ps;
      ps.allocate(4);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
      }
      // open test input file
      std::ifstream in;

      // read in data
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         in >> ps[i];
      }

      // assert that it is read correctly
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

      // set up system
      System<3> system;
      SweepTest::SetUpSystem(system, "in/diblock/bcc/param.flex");

      // set up LinearSweepParameter objects 
      DArray< LinearSweepParameter<3> > ps;
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

      // call get_ function to get value through parameter
      for (int i = 0; i < 4; ++i) {
         paramval[i] = ps[i].current();
      }

      // manually check equality for each one
      sysval[0] = system.mixture().polymer(0).block(0).length();
      sysval[1] = system.interaction().chi(0,1);
      sysval[2] = system.mixture().monomer(0).step();
      sysval[3] = system.mixture().polymer(0).phi();

      for (int i = 0; i < 4; ++i) {
         TEST_ASSERT(sysval[i]==paramval[i]);
      }
      
   }

   void testParameterSet()
   {
      printMethod(TEST_FUNC);

      // set up system
      System<3> system;
      SweepTest::SetUpSystem(system, "in/diblock/bcc/param.flex");

      // set up LinearSweepParameter objects 
      DArray< LinearSweepParameter<3> > ps;
      ps.allocate(4);
      DArray<double> initial;
      initial.allocate(4);
      std::ifstream in;
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
         in >> ps[i];
         ps[i].getInitial();
         initial[i] = ps[i].current();
      }

      DArray<double> sysval, paramval;
      sysval.allocate(4);
      paramval.allocate(4);

      // Set for some arbitrary value of s in [0,1]
      double s = 0.295586;
      for (int i = 0; i < 4; ++i) {
         ps[i].update(s);
      }
      // calculate expected value of parameter using s
      for (int i = 0; i < 4; ++i) {
         paramval[i] = initial[i] + s*ps[i].change();
      }

      // manually check equality for each one
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

      // set up system with Linear Sweep Object
      System<3> system;
      SweepTest::SetUpSystem(system, "in/param.sweep");
   }

   void testLinearSweepIterate()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepIterate");

      // Set up system with a LinearSweep object
      System<3> system;
      SweepTest::SetUpSystem(system, "in/param.sweep");
      
      // Read initial w fields and sweep
      system.readWBasis("in/diblock/bcc/omega.ref");
      system.sweep();

   }

   void SetUpSystem(System<3>& system, std::string fname)
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
//TEST_ADD(SweepTest, testLinearSweepIterate)
TEST_END(SweepTest)

#endif
