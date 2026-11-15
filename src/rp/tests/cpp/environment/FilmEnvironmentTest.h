#ifndef RPC_FILM_ENVIRONMENT_TEST_H
#define RPC_FILM_ENVIRONMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/environment/FilmEnvironment.h>
#include <rp/scft/sweep/SweepParameter.h>
#include <pscf/sweep/ParameterType.h>
#include <util/containers/GArray.h>
#include <util/math/Constants.h>
#include <util/misc/Exception.h>

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
   * Read the parameters of in/environment1 from a stream that cannot
   * seek, as for a parameter file read from a pipe.
   */
   void testReadParametersNoSeek()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestReadParametersNoSeek.log");

      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");
      Rp::FilmEnvironment<1,CPT> env(system);

      // Stream buffer with the text of in/environment1 that cannot
      // seek (std::streambuf::seekoff and seekpos return -1)
      class NoSeekBuf : public std::streambuf
      {
      public:
         NoSeekBuf(std::string const & text)
          : text_(text)
         {  setg(&text_[0], &text_[0], &text_[0] + text_.size()); }
      private:
         std::string text_;
      };
      std::ifstream file;
      openInputFile("in/environment1", file);
      std::stringstream text;
      text << file.rdbuf();
      file.close();
      NoSeekBuf buf(text.str());
      std::istream in(&buf);
      TEST_ASSERT(in.tellg() == std::streampos(-1));

      env.readParam(in);
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

   /*
   * A film without a space group, with walls that are not identical,
   * has the same mask and free energy as with space group P_1.
   */
   void testNoSpaceGroup()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestNoSpaceGroup.log");

      // Reference: the system of in/system1DEnv, with space group P_1
      Rp::System<1,CPT> ref;
      createSystem(ref, "in/system1DEnv");
      ref.w().readBasis("in/wIn1D.bf");
      ref.w().writeRGrid("out/FilmEnvTestNoSpaceGroup_w.rf");
      ref.compute();
      ref.scft().compute();

      // The same system without a space group, with AmIteratorGrid
      Rp::System<1,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } ds 0.0025 }"
         " Interaction{ chi( 1 0 20.0 ) }"
         " Domain{ mesh 96 lattice lamellar }"
         " FilmEnvironment{ normalVecId 0 interfaceThickness 0.2"
         " excludedThickness 0.4 chiBottom[ 0.0 0.0 ]"
         " chiTop[ 40.0 0.0 ] }"
         " AmIteratorGrid{ epsilon 1.0e-12 maxItr 10 maxHist 10"
         " isFlexible 0 } }\n");
      system.readParam(param);
      TEST_ASSERT(!system.domain().hasGroup());
      system.w().readRGrid("out/FilmEnvTestNoSpaceGroup_w.rf");
      system.compute();
      system.scft().compute();
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.h().hasData());
      TEST_ASSERT(eq(system.mask().phiTot(), ref.mask().phiTot()));
      double diff = std::abs(system.scft().fHelmholtz() 
                             - ref.scft().fHelmholtz());
      if (verbose() > 0) {
         std::cout << "\nFree energy difference = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-10);
   }

   /*
   * In a film with a space group but without an Iterator, the mask and
   * the wall fields are symmetric, and equal to those of the same system
   * with AmIteratorBasis, so they can be written in basis format.
   */
   void testSymmetricFieldsNoIterator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestSymmetricFieldsNoIterator.log");

      // Reference: the system of in/system1DEnv, with AmIteratorBasis
      Rp::System<1,CPT> ref;
      createSystem(ref, "in/system1DEnv");
      ref.w().readBasis("in/wIn1D.bf");
      ref.compute();
      TEST_ASSERT(ref.mask().isSymmetric());
      TEST_ASSERT(ref.h().isSymmetric());

      // The same system without an Iterator
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
      std::istringstream in(
         "READ_W_BASIS in/wIn1D.bf\nCOMPUTE\n"
         "WRITE_MASK_BASIS out/FilmEnvTestSymmetric_mask.bf\n"
         "WRITE_H_BASIS out/FilmEnvTestSymmetric_h.bf\n"
         "FINISH\n");
      system.readCommands(in);
      TEST_ASSERT(system.mask().isSymmetric());
      TEST_ASSERT(system.h().isSymmetric());

      double maxDiff = 0.0;
      int nBasis = system.domain().basis().nBasis();
      for (int k = 0; k < nBasis; ++k) {
         double diff = std::abs(system.mask().basis()[k] 
                                - ref.mask().basis()[k]);
         if (diff > maxDiff) maxDiff = diff;
         for (int i = 0; i < 2; ++i) {
            diff = std::abs(system.h().basis(i)[k] 
                            - ref.h().basis(i)[k]);
            if (diff > maxDiff) maxDiff = diff;
         }
      }
      TEST_ASSERT(maxDiff < 1.0E-12);
   }

   /*
   * A 2D film with walls normal to lattice vector a of an oblique cell.
   * A change of the angle gamma leaves vector a unchanged, but a cell
   * with gamma != pi/2 must still be rejected after a first COMPUTE at
   * gamma = pi/2, as it is when it is the first cell.
   */
   void testCellChangeCompatibility()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestCellChangeCompatibility.log");

      Rp::System<2,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } ds 0.01 }"
         " Interaction{ chi( 1 0 15.0 ) }"
         " Domain{ mesh 16 8 lattice oblique groupName p_2 }"
         " FilmEnvironment{ normalVecId 0 interfaceThickness 0.2"
         " excludedThickness 0.4 chiBottom[ 5.0 0.0 ]"
         " chiTop[ 5.0 0.0 ] } }\n");
      system.readParam(param);

      // Orthogonal cell, and w fields with only a homogeneous part
      FSArray<double, 6> parameters;
      parameters.append(3.1);
      parameters.append(1.6);
      parameters.append(0.5*Constants::Pi);
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
      w[0][0] = 1.0;
      system.w().setBasis(w);
      system.compute();

      // Change gamma, which leaves lattice vector a unchanged
      parameters[2] = 1.2;
      system.setUnitCell(parameters);
      bool thrown = false;
      try {
         system.compute();
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * WRITE_MASK_BASIS and WRITE_H_BASIS after SET_UNIT_CELL write the
   * mask and the wall fields of the new unit cell.
   */
   void testWriteFieldsAfterCellChange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestWriteFieldsAfterCellChange.log");

      // Reference: film of thickness 2.6
      Rp::System<1,CPT> ref;
      createSystem(ref, "in/system1DEnv");
      std::istringstream refCommands(
         "READ_W_BASIS in/wIn1D.bf\nSET_UNIT_CELL lamellar 2.6\n"
         "COMPUTE\nFINISH\n");
      ref.readCommands(refCommands);

      // Film of thickness 2.1 (cell of in/wIn1D.bf), then 2.6
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1DEnv");
      std::istringstream commands(
         "READ_W_BASIS in/wIn1D.bf\nCOMPUTE\nSET_UNIT_CELL lamellar 2.6\n"
         "WRITE_MASK_BASIS out/FilmEnvTestCellChange_mask.bf\n"
         "WRITE_H_BASIS out/FilmEnvTestCellChange_h.bf\nFINISH\n");
      system.readCommands(commands);

      // Compare the fields written to file with the reference
      UnitCell<1> unitCell;
      DArray<double> mask;
      DArray< DArray<double> > h;
      system.domain().fieldIo().readFieldBasis(
                       "out/FilmEnvTestCellChange_mask.bf", mask, unitCell);
      TEST_ASSERT(eq(unitCell.parameter(0), 2.6));
      system.domain().fieldIo().readFieldsBasis(
                       "out/FilmEnvTestCellChange_h.bf", h, unitCell);
      double maxDiff = 0.0;
      int nBasis = system.domain().basis().nBasis();
      for (int k = 0; k < nBasis; ++k) {
         double diff = std::abs(mask[k] - ref.mask().basis()[k]);
         if (diff > maxDiff) maxDiff = diff;
         for (int i = 0; i < 2; ++i) {
            diff = std::abs(h[i][k] - ref.h().basis(i)[k]);
            if (diff > maxDiff) maxDiff = diff;
         }
      }
      if (verbose() > 0) {
         std::cout << "\nMax difference = " << maxDiff << "\n";
      }
      TEST_ASSERT(maxDiff < 1.0E-8);
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

};

TEST_BEGIN(FilmEnvironmentTest)
TEST_ADD(FilmEnvironmentTest, testConstructor)
TEST_ADD(FilmEnvironmentTest, testReadParameters)
TEST_ADD(FilmEnvironmentTest, testReadParametersNoSeek)
TEST_ADD(FilmEnvironmentTest, testWallChiParameter)
TEST_ADD(FilmEnvironmentTest, testSolve1D)
TEST_ADD(FilmEnvironmentTest, testSolve2D)
TEST_ADD(FilmEnvironmentTest, testComputeNoIterator)
TEST_ADD(FilmEnvironmentTest, testNoSpaceGroup)
TEST_ADD(FilmEnvironmentTest, testSymmetricFieldsNoIterator)
TEST_ADD(FilmEnvironmentTest, testCellChangeCompatibility)
TEST_ADD(FilmEnvironmentTest, testWriteFieldsAfterCellChange)
TEST_ADD(FilmEnvironmentTest, testSweep)
TEST_ADD(FilmEnvironmentTest, testSolveWithFBulk)
TEST_ADD(FilmEnvironmentTest, testStressFiniteDifference1D)
TEST_END(FilmEnvironmentTest)

#endif
