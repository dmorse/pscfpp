#ifndef PSPC_SYSTEM_TEST_H
#define PSPC_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/field/BFieldComparison.h>

//#include <pscf/mesh/MeshIterator.h>
//#include <util/format/Dbl.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class SystemTest : public UnitTest
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
      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testConversion1D_lam.log");

      std::ifstream in;
      openInputFile("in/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w-fields (reference solution, solved by Fortran PSCF)
      system.readWBasis("in/diblock/lam/omega.in");

      // Copy w field components to wFields_check after reading
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in/diblock/lam/omega.in",
                          "out/testConversion1D_lam_w.rf");
      system.rGridToBasis("out/testConversion1D_lam_w.rf",
                          "out/testConversion1D_lam_w.bf");
      system.readWBasis("out/testConversion1D_lam_w.bf");

      // Compare result to original
      BFieldComparison comparison;
      comparison.compare(wFields_check, system.wFields());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

   }

   void testConversion2D_hex()
   {
      printMethod(TEST_FUNC);
      System<2> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testConversion2D_hex.log");

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in/diblock/hex/omega.in");

      // Store components in wFields_check for later comparison
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Round trip basis -> rgrid -> basis, read resulting wField
      system.basisToRGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.rf");

      system.rGridToBasis("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w.bf");
      system.readWBasis("out/testConversion2D_hex_w.bf");

      // Check symmetry of rgrid representation
      bool hasSymmetry
       = system.checkRGridFieldSymmetry("out/testConversion2D_hex_w.rf");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      BFieldComparison comparison;
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
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
      system.readWBasis("in/diblock/bcc/omega.in");

      // Store components of field as input
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Complete round trip basis -> rgrid -> basis
      system.basisToRGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.rf");
      system.rGridToBasis("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w.bf");
      system.readWBasis("out/testConversion3D_bcc_w.bf");

      // Check symmetry of rgrid representation
      bool hasSymmetry
       = system.checkRGridFieldSymmetry("out/testConversion3D_bcc_w.rf");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      BFieldComparison comparison;
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
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

      system.readWBasis("in/diblock/bcc/omega.in");
      bool hasSymmetry = system.fieldIo().hasSymmetry(system.wFieldRGrid(0));
      TEST_ASSERT(hasSymmetry);

      // Intentionally mess up the field, check that symmetry is destroyed
      system.wFieldRGrid(0)[23] += 0.1;
      hasSymmetry = system.fieldIo().hasSymmetry(system.wFieldRGrid(0));
      TEST_ASSERT(!hasSymmetry);

   }

   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_rigid.log");

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/lam/param.rigid", in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in/diblock/lam/omega.ref");

      // Make a copy of the original field
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate1D_lam_rigid_w.bf");
      system.writeCBasis("out/testIterate1D_lam_rigid_c.bf");

      // Compare solution to original
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-8) {
         stress = true;
      }
      TEST_ASSERT(stress);

   }

   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_flex.log");

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/lam/omega.ref");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate1D_lam_flex_w.bf");
      system.writeCBasis("out/testIterate1D_lam_flex_c.bf");

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_rigid.log");

      System<2> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/hex/param.rigid", in);
      system.readParam(in);
      in.close();

      // Read reference solution
      system.readWBasis("in/diblock/hex/omega.ref");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read initial guess, iterate, output solution
      system.readWBasis("in/diblock/hex/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate2D_hex_rigid_w.bf");
      system.writeCBasis("out/testIterate2D_hex_rigid_c.bf");

      // Compare current solution to reference solution
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 3.0E-7);
      // Maximum error of 2.608E-7 occurs for the first star

      // Check stress
      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-8) {
         stress = true;
      }
      TEST_ASSERT(stress);
   }

   void testIterate2D_hex_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_flex.log");

      System<2> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read reference solution (produced by Fortran code)
      system.readWBasis("in/diblock/hex/omega.ref");

      // Save reference solution to wFields_check array
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/hex/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate2D_hex_flex_w.bf");
      system.writeCBasis("out/testIterate2D_hex_flex_c.bf");

      // Compare solution to reference solution
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 3.0E-7);
      // Maximum difference of 2.58E-7 occurs for the first star

   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_rigid.log");

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/bcc/param.rigid", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/bcc/omega.ref");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate3D_bcc_rigid_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_rigid_c.bf");

      // Compare solution to reference solution
      BFieldComparison comparison(1); // Constructor argument 1 skips star 0
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum difference of 1.023E-7 occurs for the second star

      // Test that stress is small
      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-7) {
         stress = true;
      }
      TEST_ASSERT(stress);

   }

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_flex.log");

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/bcc/omega.ref");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      system.iterate();
      system.writeWBasis("out/testIterate3D_bcc_flex_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_flex_c.bf");

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum difference of 1.09288E-7 occurs for the second star

   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
TEST_ADD(SystemTest, testCheckSymmetry3D_bcc)
TEST_ADD(SystemTest, testIterate1D_lam_rigid)
TEST_ADD(SystemTest, testIterate1D_lam_flex)
TEST_ADD(SystemTest, testIterate2D_hex_rigid)
TEST_ADD(SystemTest, testIterate2D_hex_flex)
TEST_ADD(SystemTest, testIterate3D_bcc_rigid)
TEST_ADD(SystemTest, testIterate3D_bcc_flex)

TEST_END(SystemTest)

#endif
