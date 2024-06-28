#ifndef RPG_SWEEP_TEST_H
#define RPG_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/sweep/SweepFactory.h>
#include <rpg/sweep/LinearSweep.h>

#include <prdc/crystal/BFieldComparison.h>

#include <pscf/sweep/ParameterType.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/format/Dbl.h>
#include <util/containers/GArray.h>

#include <fstream>
#include <sstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class SweepTest : public LogFileUnitTest
{

public:


   void setUp()
   {  setVerbose(0); }

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
      SweepFactory<3> sf(system);
      Sweep<3>* sweepPtr;
      
      sweepPtr = sf.factory("LinearSweep");
      TEST_ASSERT(sweepPtr != 0);
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
      TEST_ASSERT(ps[0].type()=="block");
      TEST_ASSERT(ps[0].id(0)==0);
      TEST_ASSERT(ps[0].id(1)==0);
      TEST_ASSERT(ps[0].change()==0.25);
      TEST_ASSERT(ps[0].isSpecialized() == false);
      TEST_ASSERT(ps[1].type()=="chi");
      TEST_ASSERT(ps[1].id(0)==0);
      TEST_ASSERT(ps[1].id(1)==1);
      TEST_ASSERT(ps[1].change()==5.00);
      TEST_ASSERT(ps[1].isSpecialized() == false);
      TEST_ASSERT(ps[2].type()=="kuhn");
      TEST_ASSERT(ps[2].id(0)==0);
      TEST_ASSERT(ps[2].change()==0.1);
      TEST_ASSERT(ps[2].isSpecialized() == false);
      TEST_ASSERT(ps[3].type()=="phi_polymer");
      TEST_ASSERT(ps[3].id(0)==0);
      TEST_ASSERT(ps[3].change()==-0.01);
      TEST_ASSERT(ps[3].isSpecialized() == false);
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
      sysval[2] = system.mixture().monomer(0).kuhn();
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
      sysval[2] = system.mixture().monomer(0).kuhn();
      sysval[3] = system.mixture().polymer(0).phi();
      for (int i = 0; i < 4; ++i) {
         TEST_ASSERT(sysval[i]==paramval[i]);
      }
   }

   void testSpecializedParameter() // test reading of specialized parameters
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSpecializedParameter");

      // Set up ParameterType array
      GArray<ParameterType> pTypes;
      ParameterType pType1, pType2, pType3;
      ParameterModifier modifier;

      pType1.name = "test1";
      pType1.nId = 1;
      pType1.modifierPtr_ = &modifier;
      pTypes.append(pType1);

      pType2.name = "test2";
      pType2.nId = 2;
      pType2.modifierPtr_ = &modifier;
      pTypes.append(pType2);

      pType3.name = "test3";
      pType3.nId = 3;
      pType3.modifierPtr_ = &modifier;
      pTypes.append(pType3);
      Log::file() << pTypes.size() << std::endl;

      // Set up SweepParameter objects 
      DArray< SweepParameter<1> > ps;
      ps.allocate(3);
      std::ifstream in;
      openInputFile("in/special/param.test", in);
      for (int i = 0; i < 3; ++i) {
         ps[i].setParameterTypesArray(pTypes);
         in >> ps[i];
      }

      // Check that the SweepParameters were read correctly
      TEST_ASSERT(ps[0].type()=="test1");
      TEST_ASSERT(ps[0].id(0)==0);
      TEST_ASSERT(ps[0].change()==-0.08);
      TEST_ASSERT(ps[0].parameterType().name == "test1");
      TEST_ASSERT(ps[0].parameterType().nId == 1);
      TEST_ASSERT(ps[0].parameterTypeId() == 0);
      TEST_ASSERT(ps[0].isSpecialized());
      TEST_ASSERT(ps[1].type()=="test2");
      TEST_ASSERT(ps[1].id(0)==1);
      TEST_ASSERT(ps[1].id(1)==2);
      TEST_ASSERT(ps[1].change()==+0.08);
      TEST_ASSERT(ps[1].parameterType().name == "test2");
      TEST_ASSERT(ps[1].parameterType().nId == 2);
      TEST_ASSERT(ps[1].parameterTypeId() == 1);
      TEST_ASSERT(ps[1].isSpecialized());
      TEST_ASSERT(ps[2].type()=="test3");
      TEST_ASSERT(ps[2].id(0)==3);
      TEST_ASSERT(ps[2].id(1)==4);
      TEST_ASSERT(ps[2].id(2)==5);
      TEST_ASSERT(ps[2].change()==+40);
      TEST_ASSERT(ps[2].parameterType().name == "test3");
      TEST_ASSERT(ps[2].parameterType().nId == 3);
      TEST_ASSERT(ps[2].parameterTypeId() == 2);
      TEST_ASSERT(ps[2].isSpecialized());
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

      double maxDiff = testLinearSweepParam("block");
      
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepChi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepChi");

      double maxDiff = testLinearSweepParam("chi");

      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepKuhn()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepKuhn");

      double maxDiff = testLinearSweepParam("kuhn");
      
      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepPhi()   
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepPhi");

      double maxDiff = testLinearSweepParam("phi");

      TEST_ASSERT(maxDiff < 5.0e-7);
   }

   void testLinearSweepSolvent()   
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepSolvent");

      double maxDiff = testLinearSweepParam("solvent");

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

      FSArray<double, 6> parameters;
      parameters.append(1.3835);
      system.setUnitCell(parameters);
   }

   double testLinearSweepParam(std::string paramname)
   {
      // Set up system with a LinearSweep object
      System<1> system;
      SweepTest::SetUpSystem(system, "in/" + paramname + "/param");
      
      // Read expected w fields
      DArray< BasisFieldState<1> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/" + paramname + "/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.readWBasis("in/" + paramname + "/w.bf");
      system.sweep();

      // Check if sweep had to backtrack. It shouldn't need to. 
      std::ifstream f(std::string(filePrefix() + "out/" + paramname + "/5_w.bf").c_str());
      if (f.good()) {
         TEST_THROW("Sweep backtracked due to iteration count greater than maxItr.");
      }

      // Read outputted fields
      DArray< BasisFieldState<1> > fieldsOut;
      fieldsOut.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsOut[i].setSystem(system);
         fieldsOut[i].read("out/" + paramname + "/" + std::to_string(i) +"_w.bf");
      }

      // Compare output
      BFieldComparison comparison(1);
      double maxDiff = 0.0;
      for (int i = 0; i < 5; ++i) {
         comparison.compare(fieldsRef[i].fields(), fieldsOut[i].fields());
         if (comparison.maxDiff() > maxDiff) {
            maxDiff = comparison.maxDiff();
         }
      }
      //setVerbose(1);
      if (verbose() > 0) {
         Log::file() << std::endl;
         Log::file() << "maxDiff = " << Dbl(maxDiff, 14, 6) << std::endl;
      }

      return maxDiff;
   }

};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructors)
TEST_ADD(SweepTest, testFactory)
TEST_ADD(SweepTest, testParameterRead)
TEST_ADD(SweepTest, testParameterGet)
TEST_ADD(SweepTest, testParameterSet)
TEST_ADD(SweepTest, testSpecializedParameter)
TEST_ADD(SweepTest, testLinearSweepRead)
TEST_ADD(SweepTest, testLinearSweepBlock)
TEST_ADD(SweepTest, testLinearSweepChi)
TEST_ADD(SweepTest, testLinearSweepKuhn)
TEST_ADD(SweepTest, testLinearSweepPhi)
TEST_ADD(SweepTest, testLinearSweepSolvent)
TEST_END(SweepTest)

#endif
