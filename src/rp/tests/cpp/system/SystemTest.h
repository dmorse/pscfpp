#ifndef RPC_SYSTEM_TEST_H
#define RPC_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/MixtureModifier.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <rp/field/Mask.h>
#include <rp/field/FieldIo.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/crystal/BFieldComparison.h>

#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>

#include <util/misc/Exception.h>
#include <util/tests/LogFileUnitTest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* Access to the FileMaster that a ScftThermo refers to, through the
* protected member function fileMaster() of ScftThermo.
*/
template <int D>
class ScftThermoProbe : public Rp::ScftThermo<D,CPT>
{
public:

   static
   FileMaster const * fileMasterPtr(Rp::ScftThermo<D,CPT> const & thermo)
   {
      FileMaster const & (Rp::SystemConstRef<D,CPT>::*get)() const
                                    = &ScftThermoProbe::fileMaster;
      return &(thermo.*get)();
   }

};

class SystemTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {
     setVerbose(0);
     closeLogFile();
   }

   // Utility functions (used in unit tests)

   /*
   * Allocate an array of basis fields.
   */
   template <int D>
   void allocateBasisFields(Rp::System<D,CPT> const & system,
                            DArray< DArray<double> >& fields)
   {
      int nMonomer = system.mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      if (fields.isAllocated()) {
         UTIL_CHECK(fields.capacity() == nMonomer);
      } else {
         fields.allocate(nMonomer);
      }

      int nBasis = system.domain().basis().nBasis();
      UTIL_CHECK(nBasis > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (fields[i].isAllocated()) {
            UTIL_CHECK(fields[i].capacity() == nBasis);
         } else {
            fields[i].allocate(nBasis);
         }
      }

   }

   /*
   * Read basis fields into an array.
   */
   template <int D>
   void readBasisFields(Rp::System<D,CPT> const & system,
                        std::string filename,
                        DArray< DArray<double> >& fields,
                        UnitCell<D>& unitCell)
   {
      allocateBasisFields(system, fields);
      Rp::FieldIo<D,CPT> const & fieldIo = system.domain().fieldIo();
      fieldIo.readFieldsBasis(filename, fields, unitCell);
   }

   /*
   * Compare basis fields.
   */
   double compareBasis(DArray< DArray<double> > const & fields1,
                       DArray< DArray<double> > const & fields2,
                       std::string message)
   {
      BFieldComparison comparison(1);
      comparison.compare(fields1, fields2);
      double maxDiff = comparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << message << maxDiff;
      }
      return maxDiff;
   }

   /*
   * Compare system w fields to reference fields from a file.
   */
   template <int D>
   double readCompareWBasis(Rp::System<D,CPT> const & system,
                            std::string filename)
   {
      DArray< DArray<double> > fields;
      UnitCell<D> unitCell;
      readBasisFields(system, filename, fields, unitCell);
      std::string message = "max w field error = ";
      return compareBasis(fields, system.w().basis(), message);
   }

   /*
   * Compare system c fields to reference fields from a file.
   */
   template <int D>
   double readCompareCBasis(Rp::System<D,CPT> const & system,
                            std::string filename)
   {
      DArray< DArray<double> > fields;
      UnitCell<D> unitCell;
      readBasisFields(system, filename, fields, unitCell);
      std::string message = "max c field error = ";
      return compareBasis(fields, system.c().basis(), message);
   }

   /*
   * Compare Helmoltz and grand free energies to prior results
   *
   * On output, fHelhmoltz and pressure contain absolute differences
   */
   template <int D>
   void compareFreeEnergies(Rp::System<D,CPT> const & system,
                            double& fHelmholtz, double& pressure)
   {
      UTIL_CHECK(system.scft().hasData());
      fHelmholtz = std::abs(fHelmholtz - system.scft().fHelmholtz());
      pressure   = std::abs(pressure - system.scft().pressure());
      if (verbose() > 0) {
         std::cout << "\nfHelmholtz diff = " << fHelmholtz;
         std::cout << "\npressure diff   = " << pressure;
      }
   }

   /*
   * Setup a system, and read parameter file.
   */
   template <int D>
   void setupSystem(Rp::System<D,CPT>& system,
                   std::string paramFileName,
                   std::string wFileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      // Read parameter file
      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();
   }

   /*
   * Setup system - read parameter file and basis file. 
   */
   template <int D>
   void initSystem(Rp::System<D,CPT>& system,
                   std::string paramFileName,
                   std::string wFileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      // Read parameter file
      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.w().readBasis(wFileName);
   }

   /*
   * Construct file name prefix + #D_ + outSuffix
   */
   std::string makeFileRoot(std::string prefix,
                            std::string outSuffix,
                            int D)
   {
      std::string outFileRoot = prefix;
      outFileRoot += std::to_string(D) + "D_" + outSuffix;
      return outFileRoot;
   }

   /*
   * Iterate and output final fields.
   */
   template <int D>
   void iterate(Rp::System<D,CPT>& system,
                std::string const & outFileRoot)
   {
      // Iterate
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }

      // Write final fields
      system.w().writeBasis(outFileRoot + "_w.bf");
      system.c().writeBasis(outFileRoot + "_c.bf");

   }

   /*
   * Template for an iteration test, with regression testing.
   */
   template <int D>
   void testIterate(Rp::System<D,CPT>& system,
                    std::string paramFileName,
                    std::string wFileName,
                    std::string outSuffix,
                    double& wMaxDiff,
                    double& cMaxDiff,
                    bool compareCFields = true)
   {
      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterate", outSuffix, D);

      openLogFile(outFileRoot + ".log");

      initSystem(system, paramFileName, wFileName);

      iterate(system, outFileRoot);

      std::string refFileRoot;
      refFileRoot = makeFileRoot("ref/testIterate", outSuffix, D);

      wMaxDiff = readCompareWBasis(system, refFileRoot + "_w.bf");
      if (compareCFields) {
         cMaxDiff = readCompareCBasis(system, refFileRoot + "_c.bf");
      }

   }

   /*
   * Compute and return stress.
   */
   template <int D>
   FSArray<double, 6> computeStress(Rp::System<D,CPT>& system)
   {
      FSArray<double, 6> stress;
      int nParameter = system.domain().unitCell().nParameter();

      system.computeStress();
      for (int i = 0; i < nParameter; ++i) {
         stress.append(system.mixture().stress(i));
      }
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "stress =" ;
         for (int i = 0; i < nParameter; ++i) {
            std::cout << "  " << stress[i];
         }
      }
      return stress;
   }

   // Unit test functions

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Rp::System<1,CPT> system;
   }

   void testScftThermoFileMaster()
   {
      printMethod(TEST_FUNC);
      Rp::System<1,CPT> system;
      TEST_ASSERT(ScftThermoProbe<1>::fileMasterPtr(system.scft())
                  == &system.fileMaster());
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      Rp::System<1,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();
   }

   void testConversion1D_lam()
   {
      printMethod(TEST_FUNC);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testConversion","lam", 1);
      openLogFile(outFileRoot + ".log");

      Rp::System<1,CPT> system;
      initSystem(system,
                 "in/diblock/lam/param.flex", 
                 "in/diblock/lam/omega.in");
      Rp::FieldIo<1,CPT> const & fieldIo = system.domain().fieldIo();

      // Copy w field components to wFields_check after reading
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Round trip conversion basis -> rgrid -> basis, read result
      fieldIo.convertBasisToRGrid("in/diblock/lam/omega.in",
                                  "out/testConversion1D_lam_w.rf");
      fieldIo.convertRGridToBasis("out/testConversion1D_lam_w.rf",
                                  "out/testConversion1D_lam_w.bf");
      system.w().readBasis("out/testConversion1D_lam_w.bf");

      // Compare result to original
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff();
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      fieldIo.convertBasisToKGrid("in/diblock/lam/omega.in",
                                  "out/testConversion1D_lam_w.kf");
      fieldIo.convertKGridToBasis("out/testConversion1D_lam_w.kf",
                                  "out/testConversion1D_lam_w_2.bf");
      system.w().readBasis("out/testConversion1D_lam_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff();
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.w().readRGrid("out/testConversion1D_lam_w.rf");
      DArray< RField<1,CPT> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion1D_lam_w.rf",
                                  "out/testConversion1D_lam_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion1D_lam_w_2.kf",
                                  "out/testConversion1D_lam_w_2.rf");
      system.w().readRGrid("out/testConversion1D_lam_w_2.rf");

      // Compare result to original
      RFieldComparison<1,CPT> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff();
      }
      TEST_ASSERT(comparison3.maxDiff() < 1.0E-10);

   }

   /*
   * Run a command with invalid arguments, which must throw an Exception
   * without creating the output file.
   */
   template <int D>
   void checkInvalidCommand(Rp::System<D,CPT>& system,
                            std::string command,
                            std::string outFileName)
   {
      std::string path = filePrefix() + outFileName;
      std::remove(path.c_str());
      std::istringstream in(command + "\nFINISH\n");
      bool thrown = false;
      try {
         system.readCommands(in);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
      std::ifstream file(path.c_str());
      TEST_ASSERT(!file.is_open());
   }

   void testExpandReplicateArguments1D_lam()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testExpandReplicateArguments1D_lam.log");

      Rp::System<1,CPT> system;
      initSystem(system,
                 "in/diblock/lam/param.flex",
                 "in/diblock/lam/omega.in");
      system.w().writeRGrid("out/testExpandReplicate_w.rf");
      std::string in = "out/testExpandReplicate_w.rf ";

      checkInvalidCommand(system, "EXPAND_RGRID_DIMENSION " + in
                          + "out/testExpand_d4.rf 4 8 8 8",
                          "out/testExpand_d4.rf");
      checkInvalidCommand(system, "EXPAND_RGRID_DIMENSION " + in
                          + "out/testExpand_n0.rf 3 8 0",
                          "out/testExpand_n0.rf");
      checkInvalidCommand(system, "REPLICATE_UNIT_CELL " + in
                          + "out/testReplicate_0.rf 0",
                          "out/testReplicate_0.rf");
   }

   void testConversion2D_hex()
   {
      printMethod(TEST_FUNC);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testConversion","hex", 2);
      openLogFile(outFileRoot + ".log");

      Rp::System<2,CPT> system;
      initSystem(system,
                 "in/diblock/hex/param.flex", 
                 "in/diblock/hex/omega.in");
      Rp::FieldIo<2,CPT> const & fieldIo = system.domain().fieldIo();

      // Store components in wFields_check for later comparison
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Round trip basis -> rgrid -> basis, read resulting wField
      fieldIo.convertBasisToRGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.rf");

      fieldIo.convertRGridToBasis("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w.bf");
      system.w().readBasis("out/testConversion2D_hex_w.bf");

      // Check symmetry of rgrid representation
      bool hasSymmetry;
      hasSymmetry = fieldIo.hasSymmetry("out/testConversion2D_hex_w.rf");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff();
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      fieldIo.convertBasisToKGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.kf");
      fieldIo.convertKGridToBasis("out/testConversion2D_hex_w.kf",
                          "out/testConversion2D_hex_w_2.bf");
      system.w().readBasis("out/testConversion2D_hex_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff();
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.w().readRGrid("out/testConversion2D_hex_w.rf");
      DArray< RField<2,CPT> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion2D_hex_w_2.kf",
                          "out/testConversion2D_hex_w_2.rf");
      system.w().readRGrid("out/testConversion2D_hex_w_2.rf");

      // Compare result to original
      RFieldComparison<2,CPT> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff();
      }
      TEST_ASSERT(comparison3.maxDiff() < 1.0E-10);

   }

   void testConversion3D_bcc()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testConversion3D_bcc.log");

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w fields in system.wFields
      system.w().readBasis("in/diblock/bcc/omega.in");

      // Store components of field as input
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Complete round trip basis -> rgrid -> basis
      Rp::FieldIo<3,CPT> const & fieldIo = system.domain().fieldIo();
      fieldIo.convertBasisToRGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.rf");
      fieldIo.convertRGridToBasis("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w.bf");
      system.w().readBasis("out/testConversion3D_bcc_w.bf");

      // Check symmetry of rgrid representation
      bool hasSymmetry 
             = fieldIo.hasSymmetry("out/testConversion3D_bcc_w.rf");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff();
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      fieldIo.convertBasisToKGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.kf");
      fieldIo.convertKGridToBasis("out/testConversion3D_bcc_w.kf",
                          "out/testConversion3D_bcc_w_2.bf");
      system.w().readBasis("out/testConversion3D_bcc_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff();
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.w().readRGrid("out/testConversion3D_bcc_w.rf");
      DArray< RField<3,CPT> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion3D_bcc_w_2.kf",
                          "out/testConversion3D_bcc_w_2.rf");
      system.w().readRGrid("out/testConversion3D_bcc_w_2.rf");

      // Compare result to original
      RFieldComparison<3,CPT> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff();
      }
      TEST_ASSERT(comparison3.maxDiff() < 1.0E-10);

   }

   void testCheckSymmetry3D_bcc()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testSymmetry3D_bcc.log");

      // Read system parameter file
      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      system.w().readBasis("in/diblock/bcc/omega.in");
      bool hasSymmetry
            = system.domain().fieldIo().hasSymmetry(system.w().rgrid(0));
      TEST_ASSERT(hasSymmetry);

      // Copy the wFieldsRGrid to a temporary container
      RField<3,CPT> field;
      field.allocate(system.domain().mesh().dimensions());
      int meshSize = system.domain().mesh().size();
      for (int j = 0; j < meshSize; ++j) {
         field[j] = system.w().rgrid(0)[j];
      }
      hasSymmetry = system.domain().fieldIo().hasSymmetry(field);
      TEST_ASSERT(hasSymmetry);

      // Intentionally mess up the field, check that symmetry is destroyed
      field[23] += 0.1;
      hasSymmetry = system.domain().fieldIo().hasSymmetry(field);
      TEST_ASSERT(!hasSymmetry);

   }

   /*
   * Command scripts that end without a newline, or with an argument
   * that cannot be read.
   */
   void testReadCommandsInput()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReadCommandsInput.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/diblock/lam/param.flex", "");

      // A last command without a newline must be executed
      std::istringstream in1("READ_W_BASIS in/diblock/lam/omega.in\n"
                             "COMPUTE");
      system.readCommands(in1);
      TEST_ASSERT(system.c().hasData());

      // A missing argument at the end of the script must be reported
      std::istringstream in2("WRITE_Q_TAIL out/testReadCommandsInput_q 0 0");
      bool thrown = false;
      try {
         system.readCommands(in2);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      // So must an argument that cannot be read
      std::istringstream in3("WRITE_Q_TAIL out/testReadCommandsInput_q"
                             " 0 x 1\nFINISH\n");
      thrown = false;
      try {
         system.readCommands(in3);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * Commands that write solver output after a change of the w fields
   * or unit cell, without a new solution.
   */
   void testWriteNeedsSolution()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testWriteNeedsSolution.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/diblock/lam/param.flex", "");
      std::istringstream in0("READ_W_BASIS in/diblock/lam/omega.in\n"
                             "COMPUTE\nFINISH\n");
      system.readCommands(in0);

      std::string out = " out/testWriteNeedsSolution";
      std::string commands[] = {"WRITE_STRESS" + out + "_stress",
                                "WRITE_BLOCK_C_RGRID" + out + "_c.rf",
                                "WRITE_Q_SLICE" + out + "_q 0 0 0 0",
                                "WRITE_Q_TAIL" + out + "_q 0 0 0",
                                "WRITE_Q" + out + "_q 0 0 0",
                                "WRITE_Q_ALL" + out + "_qall"};
      std::string changes[] = {"READ_W_BASIS in/diblock/lam/omega.in",
                               "SET_UNIT_CELL lamellar 1.4"};
      for (int i = 0; i < 6; ++i) {

         // With a solution for the current w fields and unit cell
         std::istringstream in1("COMPUTE\n" + commands[i] + "\nFINISH\n");
         system.readCommands(in1);

         // After a change, without a new solution
         for (int j = 0; j < 2; ++j) {
            std::istringstream in2(changes[j] + "\n" + commands[i]
                                   + "\nFINISH\n");
            bool thrown = false;
            try {
               system.readCommands(in2);
            } catch (Exception& e) {
               thrown = true;
            }
            TEST_ASSERT(thrown);
         }
      }
   }

   /*
   * A parameter file written by WRITE_PARAM for a bead model system
   * must be read as a bead model system.
   */
   void testWriteParamBead()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testWriteParamBead.log");

      {
         Rp::System<1,CPT> system;
         setupSystem(system, "in/diblock/lam_bead/param.flex", "");
         TEST_ASSERT(PolymerModel::isBead());
         std::istringstream in("WRITE_PARAM out/testWriteParamBead.prm\n"
                               "FINISH\n");
         system.readCommands(in);
      }

      // The file must begin with "System{ polymerModel bead"
      std::ifstream file;
      openInputFile("out/testWriteParamBead.prm", file);
      std::string label, value;
      file >> label >> label >> value;
      file.close();
      TEST_ASSERT(label == "polymerModel");
      TEST_ASSERT(value == "bead");

      Rp::System<1,CPT> system;
      setupSystem(system, "out/testWriteParamBead.prm", "");
      TEST_ASSERT(PolymerModel::isBead());
      TEST_ASSERT(system.mixture().polymer(0).nBead() == 100);
   }

   /*
   * Status returned by readCommands: 0 if no command failed, 1 after
   * an unknown command or an ITERATE command that did not converge.
   */
   void testReadCommandsStatus()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReadCommandsStatus.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/diblock/lam/param.flex", "");
      std::istringstream in1("READ_W_BASIS in/diblock/lam/omega.in\n"
                             "COMPUTE\nFINISH\n");
      TEST_ASSERT(system.readCommands(in1) == 0);
      std::istringstream in2("COMPUTE\n");
      TEST_ASSERT(system.readCommands(in2) == 0);
      std::istringstream in3("COMPUTEE\nFINISH\n");
      TEST_ASSERT(system.readCommands(in3) == 1);

      // An AmIterator with maxItr 1 cannot converge
      Rp::System<1,CPT> system2;
      system2.fileMaster().setInputPrefix(filePrefix());
      system2.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } ds 0.01 }"
         " Interaction{ chi( 1 0 15.0 ) }"
         " Domain{ mesh 32 lattice lamellar groupName P_-1 }"
         " AmIteratorBasis{ epsilon 1.0e-10 maxItr 1 maxHist 10 } }\n");
      system2.readParam(param);
      std::istringstream in4("READ_W_BASIS in/diblock/lam/omega.in\n"
                             "ITERATE\nFINISH\n");
      TEST_ASSERT(system2.readCommands(in4) == 1);
   }

   /*
   * ANALYZE in a system without a Simulator must throw an Exception.
   */
   void testAnalyzeNoSimulator()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeNoSimulator.log");

      Rp::System<1,CPT> system;
      setupSystem(system, "in/diblock/lam/param.flex", "");
      TEST_ASSERT(!system.hasSimulator());
      std::istringstream in("ANALYZE 0 1 RGridTrajectoryReader"
                            " in/diblock/lam/omega.in\nFINISH\n");
      bool thrown = false;
      try {
         system.readCommands(in);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * Free energy and pressure of homogeneous blends of two homopolymers
   * of length 1 with volume fractions 1 - phi and phi, in fields w = 0
   * (basis) and w = 20 (r-grid), which have exact values. The terms of
   * phi = 5.0E-9 and of chi = 1.0E-9 must be included.
   */
   void testThermoSmallTerms()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testThermoSmallTerms.log");

      double phis[] = {5.0E-9, 0.5};
      double chis[] = {10.0, 1.0E-9};
      for (int k = 0; k < 2; ++k) {
         double phi = phis[k];
         double chi = chis[k];
         std::ostringstream text;
         text << std::setprecision(17) << "System{ Mixture{ nMonomer 2"
              << " monomers[ 1.0 1.0 ] nPolymer 2"
              << " Polymer{ type linear nBlock 1 blocks[ 0 1.0 ] phi "
              << 1.0 - phi << " } Polymer{ type linear nBlock 1"
              << " blocks[ 1 1.0 ] phi " << phi << " } ds 0.05 }"
              << " Interaction{ chi( 1 0 " << chi << " ) }"
              << " Domain{ mesh 8 lattice lamellar groupName P_-1 } }\n";
         Rp::System<1,CPT> system;
         std::istringstream in(text.str());
         system.readParam(in);
         FSArray<double, 6> parameters;
         parameters.append(1.0);
         system.setUnitCell(parameters);
         double fExact = (1.0 - phi)*(std::log(1.0 - phi) - 1.0)
                       + phi*(std::log(phi) - 1.0) + chi*phi*(1.0 - phi);
         DArray< RField<1,CPT> > w;
         w.allocate(2);
         w[0].allocate(system.domain().mesh().dimensions());
         w[1].allocate(system.domain().mesh().dimensions());
         for (int j = 0; j < 2; ++j) {
            for (int m = 0; m < 8; ++m) {
               w[0][m] = w[1][m] = 20.0*j;
            }
            // Symmetric (basis) fields for w = 0, r-grid fields for 20
            system.w().setRGrid(w, j == 0);
            system.compute();
            system.scft().compute();
            double f = system.scft().fHelmholtz();
            double p = system.scft().pressure();
            double pExact = 20.0*j + 1.0 - chi*phi*(1.0 - phi);
            TEST_ASSERT(std::abs(f - fExact) < 2.0E-11);
            TEST_ASSERT(std::abs(p - pExact) < 2.0E-11);
         }
      }
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testScftThermoFileMaster)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testExpandReplicateArguments1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
TEST_ADD(SystemTest, testCheckSymmetry3D_bcc)
TEST_ADD(SystemTest, testThermoSmallTerms)
TEST_ADD(SystemTest, testReadCommandsInput)
TEST_ADD(SystemTest, testWriteNeedsSolution)
TEST_ADD(SystemTest, testWriteParamBead)
TEST_ADD(SystemTest, testReadCommandsStatus)
TEST_ADD(SystemTest, testAnalyzeNoSimulator)
TEST_END(SystemTest)

#endif
