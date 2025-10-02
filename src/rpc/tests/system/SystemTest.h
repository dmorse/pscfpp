#ifndef RPC_SYSTEM_TEST_H
#define RPC_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/system/System.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>
#include <util/format/Dbl.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

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
   void allocateBasisFields(System<D> const & system,
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
   void readBasisFields(System<D> const & system,
                        std::string filename,
                        DArray< DArray<double> >& fields,
                        UnitCell<D>& unitCell)
   {
      allocateBasisFields(system, fields);
      FieldIo<D> const & fieldIo = system.domain().fieldIo();
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
   double readCompareWBasis(System<D> const & system,
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
   double readCompareCBasis(System<D> const & system,
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
   void compareFreeEnergies(System<D> const & system,
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
   void setupSystem(System<D>& system,
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
   void initSystem(System<D>& system,
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
   void iterate(System<D>& system,
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
   void testIterate(System<D>& system,
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
   FSArray<double, 6> computeStress(System<D>& system)
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
      System<1> system;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;
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

      System<1> system;
      initSystem(system,
                 "in/diblock/lam/param.flex", 
                 "in/diblock/lam/omega.in");
      FieldIo<1> const & fieldIo = system.domain().fieldIo();

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
      DArray< RField<1> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion1D_lam_w.rf",
                                  "out/testConversion1D_lam_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion1D_lam_w_2.kf",
                                  "out/testConversion1D_lam_w_2.rf");
      system.w().readRGrid("out/testConversion1D_lam_w_2.rf");

      // Compare result to original
      RFieldComparison<1> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff();
      }
      TEST_ASSERT(comparison3.maxDiff() < 1.0E-10);

   }

   void testConversion2D_hex()
   {
      printMethod(TEST_FUNC);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testConversion","hex", 2);
      openLogFile(outFileRoot + ".log");

      System<2> system;
      initSystem(system,
                 "in/diblock/hex/param.flex", 
                 "in/diblock/hex/omega.in");
      FieldIo<2> const & fieldIo = system.domain().fieldIo();

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
      DArray< RField<2> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion2D_hex_w_2.kf",
                          "out/testConversion2D_hex_w_2.rf");
      system.w().readRGrid("out/testConversion2D_hex_w_2.rf");

      // Compare result to original
      RFieldComparison<2> comparison3;
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
      System<3> system;
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
      FieldIo<3> const & fieldIo = system.domain().fieldIo();
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
      DArray< RField<3> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      fieldIo.convertRGridToKGrid("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w_2.kf");
      fieldIo.convertKGridToRGrid("out/testConversion3D_bcc_w_2.kf",
                          "out/testConversion3D_bcc_w_2.rf");
      system.w().readRGrid("out/testConversion3D_bcc_w_2.rf");

      // Compare result to original
      RFieldComparison<3> comparison3;
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
      System<3> system;
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
      RField<3> field;
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

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
TEST_ADD(SystemTest, testCheckSymmetry3D_bcc)
TEST_END(SystemTest)

#endif
