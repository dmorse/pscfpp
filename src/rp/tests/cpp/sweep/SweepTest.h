#ifndef RPC_SWEEP_TEST_H
#define RPC_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/sweep/SweepFactory.h>
#include <rp/scft/sweep/LinearSweep.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <pscf/crystal/BFieldComparison.h>

#include <pscf/sweep/ParameterType.h>
#include <pscf/interaction/Interaction.h>

#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <util/containers/GArray.h>
#include <util/format/Dbl.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <vector>

using namespace Util;
using namespace Pscf;

/*
* LinearSweep for tests of the sequence of states of a sweep.
*
* Records the value of s of each attempted state, and replaces the
* iterative solution of each state by a solution of the MDE for the
* current w fields. Solve number failSolve (counting from 0 for the
* initial state) is reported as a failure to converge.
*/
template <int D>
class TestLinearSweep : public Rp::LinearSweep<D,CPT>
{

public:

   TestLinearSweep(Rp::System<D,CPT>& system)
    : Rp::LinearSweep<D,CPT>(system),
      failSolve(-1),
      nSolve_(0)
   {}

   void setParameters(double s) override
   {
      sAttempt.push_back(s);
      Rp::LinearSweep<D,CPT>::setParameters(s);
   }

   /// Values of s of all attempted states after the initial one.
   std::vector<double> sAttempt;

   /// Index of the solve to report as failure (-1 for none).
   int failSolve;

protected:

   int solve(bool isContinuation) override
   {
      this->system().compute();
      int error = (nSolve_ == failSolve) ? 1 : 0;
      ++nSolve_;
      return error;
   }

private:

   /// Number of calls to solve.
   int nSolve_;

};

class SweepTest : public LogFileUnitTest
{

public:


   void setUp()
   {  setVerbose(0); }

   void testConstructors()
   {
      printMethod(TEST_FUNC);

      Rp::System<3,CPT> system;
      Rp::LinearSweep<3,CPT> ls(system);
      Rp::SweepFactory<3,CPT> sf(system);
   }

   void testFactory() 
   {
      printMethod(TEST_FUNC);

      Rp::System<3,CPT> system;
      Rp::SweepFactory<3,CPT> sf(system);
      Rp::Sweep<3,CPT>* sweepPtr;
      
      sweepPtr = sf.factory("LinearSweep");
      TEST_ASSERT(sweepPtr != 0);
   }

   void testParameterRead() 
   {
      printMethod(TEST_FUNC);

      // Set up system with some data
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");
      // Set up SweepParameter objects 
      DArray< Rp::SweepParameter<1,CPT> > ps;
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
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");

      // Set up SweepParameter objects 
      DArray< Rp::SweepParameter<1,CPT> > ps;
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
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");

      // Set up SweepParameter objects 
      DArray< Rp::SweepParameter<1,CPT> > ps;
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
      pType1.modifierPtr = &modifier;
      pTypes.append(pType1);

      pType2.name = "test2";
      pType2.nId = 2;
      pType2.modifierPtr = &modifier;
      pTypes.append(pType2);

      pType3.name = "test3";
      pType3.nId = 3;
      pType3.modifierPtr = &modifier;
      pTypes.append(pType3);
      Log::file() << pTypes.size() << std::endl;

      // Set up SweepParameter objects 
      DArray< Rp::SweepParameter<1,CPT> > ps;
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

   void testParameterIdRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testParameterIdRange");

      // Diblock melt, 1 lattice parameter. Parameters with an index out
      // of range must be rejected by getInitial, with a message that
      // says so.
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");
      std::string bad[] = {"phi_polymer 1 0.1", "cell_param 1 0.1",
                           "block 0 2 0.1"};
      for (int i = 0; i < 3; ++i) {
         Rp::SweepParameter<1,CPT> parameter;
         parameter.setSystem(system);
         std::istringstream in(bad[i]);
         in >> parameter;
         std::string message;
         try {
            parameter.getInitial();
         } catch (Exception& e) {
            message = e.message();
         }
         TEST_ASSERT(message.find("out of range") != std::string::npos);
      }
   }

   void testLinearSweepRead()
   {
      printMethod(TEST_FUNC);

      // Set up system with Linear Sweep Object
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");
   }

   void testLinearSweepReadNs()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepReadNs");

      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");

      // The number of steps ns and historyCapacity must be positive
      std::string bad[] = {"ns 0", "ns -2", "ns 4 historyCapacity 0"};
      for (int i = 0; i < 3; ++i) {
         Rp::LinearSweep<1,CPT> sweep(system);
         std::istringstream in(bad[i]
                               + " nParameter 1 parameters chi 0 1 1.0\n");
         bool thrown = false;
         try {
            sweep.readParameters(in);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown);
      }

      // Positive values
      Rp::LinearSweep<1,CPT> sweep(system);
      std::istringstream in("ns 4 historyCapacity 2"
                            " nParameter 1 parameters chi 0 1 1.0\n");
      sweep.readParameters(in);
   }

   void testLinearSweepPhiMixed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepPhiMixed");

      // Diblock (closed) in an open solvent: a sweep of the volume
      // fraction of the diblock alone must be accepted
      Rp::System<1,CPT> system;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1 nSolvent 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 0.8 }"
         " Solvent{ monomerId 1 size 0.1 mu 1.0 } ds 0.01 }"
         " Interaction{ chi( 1 0 12.0 ) }"
         " Domain{ mesh 40 lattice lamellar groupName P_-1 }"
         " AmIteratorBasis{ epsilon 1.0e-8 maxItr 100 maxHist 10 } }\n");
      system.readParam(param);
      Rp::LinearSweep<1,CPT> sweep(system);
      std::istringstream in("ns 2 nParameter 1"
                            " parameters phi_polymer 0 -0.1\n");
      bool thrown = false;
      try {
         sweep.readParameters(in);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
   }

   void testLinearSweepExactS()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearSweepExactS");

      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/block/param");
      system.w().readBasis("in/block/w.bf");

      // Without failures, the states must be at s = k/ns exactly, and
      // the last state at s = 1 (with ns = 9, repeated addition of 1/9
      // gives 1 + 2.2E-16)
      {
         TestLinearSweep<1> sweep(system);
         std::istringstream in("ns 9"
                               " baseFileName out/testLinearSweepExactS_"
                               " nParameter 1 parameters chi 0 1 0.9\n");
         sweep.readParameters(in);
         sweep.sweep();
         TEST_ASSERT(sweep.sAttempt.size() == 9);
         for (int k = 1; k <= 9; ++k) {
            TEST_ASSERT(sweep.sAttempt[k-1] == double(k)/9.0);
         }
      }

      // With ns = 3 and a failure at s = 2/3, the step is halved, and
      // the states must be at exact multiples of 1/6, the last at s = 1
      {
         system.w().readBasis("in/block/w.bf");
         TestLinearSweep<1> sweep(system);
         sweep.failSolve = 2;
         std::istringstream in("ns 3"
                               " baseFileName out/testLinearSweepExactS_"
                               " nParameter 1 parameters chi 0 1 0.9\n");
         sweep.readParameters(in);
         sweep.sweep();
         double s[] = {1.0/3.0, 2.0/3.0, 3.0/6.0, 4.0/6.0, 5.0/6.0, 1.0};
         TEST_ASSERT(sweep.sAttempt.size() == 6);
         for (int k = 0; k < 6; ++k) {
            TEST_ASSERT(sweep.sAttempt[k] == s[k]);
         }
      }
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

   void SetUpSystem(Rp::System<1,CPT>& system, std::string fname)
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
      Rp::System<1,CPT> system;
      SweepTest::SetUpSystem(system, "in/" + paramname + "/param");
      
      // Read expected w fields
      DArray< Rp::BasisFieldState<1,CPT> > fieldsRef;
      fieldsRef.allocate(5);
      for (int i = 0; i < 5; ++i) {
         fieldsRef[i].setSystem(system);
         fieldsRef[i].read("in/sweepref/" + paramname + "/" + std::to_string(i) +"_w.bf");
      }

      // Read initial field guess and sweep
      system.w().readBasis("in/" + paramname + "/w.bf");
      system.sweep();

      // Check if sweep had to backtrack. It shouldn't need to. 
      std::ifstream f(std::string(filePrefix() + "out/" + paramname + "/5_w.bf").c_str());
      if (f.good()) {
         TEST_THROW("Sweep backtracked due to iteration count greater than maxItr.");
      }

      // Read outputted fields
      DArray< Rp::BasisFieldState<1,CPT> > fieldsOut;
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
TEST_ADD(SweepTest, testParameterIdRange)
TEST_ADD(SweepTest, testLinearSweepRead)
TEST_ADD(SweepTest, testLinearSweepReadNs)
TEST_ADD(SweepTest, testLinearSweepPhiMixed)
TEST_ADD(SweepTest, testLinearSweepExactS)
TEST_ADD(SweepTest, testLinearSweepBlock)
TEST_ADD(SweepTest, testLinearSweepChi)
TEST_ADD(SweepTest, testLinearSweepKuhn)
TEST_ADD(SweepTest, testLinearSweepPhi)
TEST_ADD(SweepTest, testLinearSweepSolvent)
TEST_END(SweepTest)

#endif
