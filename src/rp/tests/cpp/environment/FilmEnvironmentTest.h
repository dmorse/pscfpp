#ifndef RPC_FILM_ENVIRONMENT_TEST_H
#define RPC_FILM_ENVIRONMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/environment/FilmEnvironment.h>
#include <rp/scft/sweep/SweepParameter.h>
#include <pscf/sweep/ParameterType.h>
#include <util/containers/GArray.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/FieldGenerator.h>
#include <util/format/Dbl.h>
#include <util/containers/FSArray.h>
#include <rp/solvers/Mixture.h>
#include <rp/environment/MixAndMatchEnv.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <pscf/chem/Ensemble.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/BFieldComparison.h>

#include <sstream>
#include <fstream>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Rp;

class FilmEnvironmentTest : public UnitTest
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
      setVerbose(0);
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Rp::System<1,CPT> system;
      Rp::FilmEnvironment<1,CPT> ext(system);
   }

   void testReadParameters() // test FilmEnvironment::readParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestReadParameters.log");

      // Set up film environment from file
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");
      Rp::FilmEnvironment<1,CPT> env(system);

      std::ifstream in;
      openInputFile("in/environment1", in);
      env.readParam(in);
      in.close();

      // Check that the parameters that are publicly accessible were 
      // read correctly
      TEST_ASSERT(env.fieldGenerator1().type() == FieldGenerator::Mask);
      TEST_ASSERT(env.fieldGenerator2().type() == FieldGenerator::External);
      
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 0;
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_bottom",ids), 5.0));
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_top",ids), 2.0));
      ids[0] = 1;
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_bottom",ids), 0.0));
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_top",ids), 10.0));
   }

   /*
   * A change of chi_bottom through the sweep parameter types of the
   * Environment, as made by a LinearSweep, must reset the Environment,
   * so that the wall fields are regenerated.
   */
   void testWallChiParameter()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestWallChiParameter.log");

      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");
      system.w().readBasis("in/wIn1D.bf");
      system.compute();
      TEST_ASSERT(!system.environment().needsUpdate());

      GArray<ParameterType> types = system.environment().getParameterTypes();
      Rp::SweepParameter<1,CPT> parameter;
      parameter.setSystem(system);
      parameter.setParameterTypesArray(types);
      std::istringstream in("chi_bottom 0 5.0");
      in >> parameter;
      parameter.getInitial();
      parameter.update(parameter.initial() + 5.0);
      TEST_ASSERT(system.environment().needsUpdate());
   }

   void testSolve1D() // solve a 1D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve1D.log");
      
      // Set up system with some data
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef1D.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      system.domain().fieldIo().writeFieldsBasis("out/w1D.bf", 
                                                 system.w().basis(), 
                                                 system.domain().unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-2);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.1117881919) < 1e-4);

      // output mask field for reference
      system.domain().fieldIo().writeFieldBasis("out/mask_1D.bf", 
                                                system.mask().basis(),
                                                system.domain().unitCell());
      
      // output external field for reference
      system.domain().fieldIo().writeFieldsBasis("out/h_1D.bf", 
                                                 system.h().basis(),
                                                 system.domain().unitCell());
   }

   void testSolve2D() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2D.log");
      
      // Set up system with some data
      Rp::System<2,CPT> system;
      createSystem(system, "in/system2DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      system.w().writeBasis("out/w2D.bf");
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536665618;
      TEST_ASSERT(std::abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef2D.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
     
      double epsilon = 1.0E-2; 
      double diff = bComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.91037539493) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.8397347928) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.91037539493) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.8397347928) < 1e-4);
   }

   /*
   * A film system without an Iterator: COMPUTE must generate the mask
   * and the wall fields and solve the MDE, and WRITE_STRESS must work.
   */
   void testComputeNoIterator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestComputeNoIterator.log");

      // The system of in/system1DEnv, without an Iterator
      Rp::System<1,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } ds 0.0025 }"
         " Interaction{ chi( 1 0 20.0 ) }"
         " Domain{ mesh 96 lattice lamellar groupName P_1 }"
         " FilmEnvironment{ normalVecId 0 interfaceThickness 0.2"
         " excludedThickness 0.4 chiBottom[ 0.0 0.0 ]"
         " chiTop[ 40.0 0.0 ] } }\n");
      system.readParam(param);
      TEST_ASSERT(!system.hasIterator());

      std::istringstream in("READ_W_BASIS in/wIn1D.bf\nCOMPUTE\n"
                            "WRITE_STRESS out/FilmEnvTestComputeNoIterator\n"
                            "FINISH\n");
      system.readCommands(in);
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(system.c().hasData());
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweep.log");
      
      // Set up system
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("out/w1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRefSweep.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      double diff = bComparison.maxDiff();

      double epsilon = 1.0E-2; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.87318676998) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.0498211637) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulk.log");
      
      // Set up system with some data
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnvFBulk");

      // Read initial guess
      system.w().readBasis("in/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramVal = system.domain().unitCell().parameter(0);
      double paramRef = 2.06125919708e+00; // changed for #272
      double paramErr = std::abs(paramVal - paramRef);
      //double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      if (verbose() > 0) {
         std::cout << "\nCell parameter value = " << Dbl(paramVal, 12, 20);
         std::cout << "\nCell parameter diff  = " << Dbl(paramErr, 12, 20);
      }
      TEST_ASSERT(paramErr < 1.0e-5);
      TEST_ASSERT(std::abs(system.mask().phiTot() - 0.8059299672) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef1D_2.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      system.domain().fieldIo().writeFieldsBasis("out/w1D_fBulk.bf", 
                                                 system.w().basis(), 
                                                 system.domain().unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.80031009905) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.9405443827) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.80031009905) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.9405443827) < 1e-4);
   }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(Rp::System<D,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }


   /*
   * Compare the stress of a 1D film (bead model, rigid cell, athermal
   * walls) with the derivative of fHelmholtz with respect to the film
   * thickness L, computed by central finite differences (issue #272).
   * The stress is the sum of the Mixture stress and the stress terms
   * of the field generators. The mask term contains d(phiTot)/dL.
   */
   void testStressFiniteDifference1D()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestStressFiniteDifference1D.log");

      // Converged solution for a film of thickness L
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DStressFD");
      system.w().readBasis("in/w1DStressFD.bf");
      double L = system.domain().unitCell().parameter(0);
      int error = system.iterate();
      TEST_ASSERT(!error);

      // Stress, as used by the iterator for a flexible film thickness
      system.computeStress();
      Rp::MixAndMatchEnv const & env
            = dynamic_cast<Rp::MixAndMatchEnv const &>(system.environment());
      double stress = system.mixture().stress(0);
      stress += env.fieldGenerator1().stress(0);
      if (env.hasFieldGenerator2()) {
         stress += env.fieldGenerator2().stress(0);
      }

      // Free energies for thicknesses L + h and L - h
      double h = 1.0E-4;
      FSArray<double, 6> parameters;
      parameters.append(L + h);
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double fp = system.scft().fHelmholtz();

      system.w().readBasis("in/w1DStressFD.bf");
      parameters[0] = L - h;
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double fm = system.scft().fHelmholtz();

      double stressFD = (fp - fm)/(2.0*h);
      double diff = std::abs(stress - stressFD);
      Log::file() << "\nstress            = " << Dbl(stress, 20, 12)
                  << "\nstress (fin diff) = " << Dbl(stressFD, 20, 12)
                  << "\ndifference        = " << Dbl(diff, 20, 12)
                  << std::endl;
      TEST_ASSERT(diff < 1.0E-7);
   }


   /*
   * Semi-grand potential g = fHelmholtz - sum of phi*mu/N over the open
   * species, for the bead model.
   */
   double semiGrandPotential(Rp::System<1,CPT> const & system)
   {
      Rp::Mixture<1,CPT> const & mixture = system.mixture();
      double g = system.scft().fHelmholtz();
      for (int i = 0; i < mixture.nPolymer(); ++i) {
         if (mixture.polymer(i).ensemble() == Ensemble::Open) {
            g -= mixture.polymer(i).phi() * mixture.polymer(i).mu()
                 / double(mixture.polymer(i).nBead());
         }
      }
      for (int i = 0; i < mixture.nSolvent(); ++i) {
         if (mixture.solvent(i).ensemble() == Ensemble::Open) {
            g -= mixture.solvent(i).phi() * mixture.solvent(i).mu()
                 / mixture.solvent(i).size();
         }
      }
      return g;
   }

   /*
   * A 1D film (bead model, fBulk given) with a closed diblock and an
   * open homopolymer and solvent. The stress modified by the mask, which
   * is used to optimize the film thickness L, must be the derivative of
   * (L - t)*(g - fBulk) with respect to L, where t is the excluded
   * thickness and g is the semi-grand potential, as computed by central
   * finite differences of solutions for a rigid unit cell.
   */
   void testModifyStressOpenSpecies()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestModifyStressOpenSpecies.log");

      Rp::System<1,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ polymerModel bead"
         " Mixture{ nMonomer 2 monomers[ 0.2236067977 0.2236067977 ]"
         " nPolymer 2 nSolvent 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 8 1 12 ] phi 0.8 }"
         " Polymer{ type linear nBlock 1 blocks[ 0 10 ] mu 0.27 }"
         " Solvent{ monomerId 1 size 2.0 mu -1.7 } }"
         " Interaction{ chi( 1 0 1.0 ) }"
         " Domain{ mesh 48 lattice lamellar groupName P_-1 }"
         " FilmEnvironment{ normalVecId 0 interfaceThickness 0.2"
         " excludedThickness 0.4 fBulk 0.03"
         " chiBottom[ 0.5 0.0 ] chiTop[ 0.5 0.0 ] }"
         " AmIteratorBasis{ epsilon 1.0e-12 maxItr 1000 maxHist 25"
         " isFlexible 0 } }\n");
      system.readParam(param);
      double t = 0.4;
      double fBulk = 0.03;

      // Solution for L = 2.2, from homogeneous w fields
      double L = 2.2;
      FSArray<double, 6> parameters;
      parameters.append(L);
      system.setUnitCell(parameters);
      int nBasis = system.domain().basis().nBasis();
      DArray< DArray<double> > w;
      w.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(nBasis);
         for (int k = 0; k < nBasis; ++k) {
            w[i][k] = 0.0;
         }
      }
      w[0][0] = 0.6;
      w[1][0] = 0.4;
      system.w().setBasis(w);
      int error = system.iterate();
      TEST_ASSERT(!error);
      w = system.w().basis();

      // Stress modified by the mask
      system.computeStress();
      Rp::MixAndMatchEnv const & env
            = dynamic_cast<Rp::MixAndMatchEnv const &>(system.environment());
      double stress = system.mixture().stress(0);
      stress += env.fieldGenerator1().stress(0);
      stress += env.fieldGenerator2().stress(0);
      stress = env.fieldGenerator1().modifyStress(0, stress);

      // (L - t)*(g - fBulk) for thicknesses L + h and L - h
      double h = 1.0E-4;
      parameters[0] = L + h;
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double gp = (L + h - t) * (semiGrandPotential(system) - fBulk);

      system.w().setBasis(w);
      parameters[0] = L - h;
      system.setUnitCell(parameters);
      error = system.iterate();
      TEST_ASSERT(!error);
      double gm = (L - h - t) * (semiGrandPotential(system) - fBulk);

      double stressFD = (gp - gm)/(2.0*h);
      double diff = std::abs(stress - stressFD);
      Log::file() << "\nmodified stress   = " << Dbl(stress, 20, 12)
                  << "\nstress (fin diff) = " << Dbl(stressFD, 20, 12)
                  << "\ndifference        = " << Dbl(diff, 20, 12)
                  << std::endl;
      TEST_ASSERT(diff < 1.0E-6);
   }

};

TEST_BEGIN(FilmEnvironmentTest)
TEST_ADD(FilmEnvironmentTest, testConstructor)
TEST_ADD(FilmEnvironmentTest, testReadParameters)
TEST_ADD(FilmEnvironmentTest, testWallChiParameter)
TEST_ADD(FilmEnvironmentTest, testSolve1D)
TEST_ADD(FilmEnvironmentTest, testSolve2D)
TEST_ADD(FilmEnvironmentTest, testComputeNoIterator)
TEST_ADD(FilmEnvironmentTest, testSweep)
TEST_ADD(FilmEnvironmentTest, testSolveWithFBulk)
TEST_ADD(FilmEnvironmentTest, testStressFiniteDifference1D)
TEST_ADD(FilmEnvironmentTest, testModifyStressOpenSpecies)
TEST_END(FilmEnvironmentTest)

#endif
