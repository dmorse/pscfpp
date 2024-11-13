#ifndef RPC_SYSTEM_TEST_H
#define RPC_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>

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
      wFields_check = system.w().basis();

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in/diblock/lam/omega.in",
                          "out/testConversion1D_lam_w.rf");
      system.rGridToBasis("out/testConversion1D_lam_w.rf",
                          "out/testConversion1D_lam_w.bf");
      system.readWBasis("out/testConversion1D_lam_w.bf");

      // Compare result to original
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      system.basisToKGrid("in/diblock/lam/omega.in",
                          "out/testConversion1D_lam_w.kf");
      system.kGridToBasis("out/testConversion1D_lam_w.kf",
                          "out/testConversion1D_lam_w_2.bf");
      system.readWBasis("out/testConversion1D_lam_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.readWRGrid("out/testConversion1D_lam_w.rf");
      DArray< RField<1> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      system.rGridToKGrid("out/testConversion1D_lam_w.rf",
                          "out/testConversion1D_lam_w_2.kf");
      system.kGridToRGrid("out/testConversion1D_lam_w_2.kf",
                          "out/testConversion1D_lam_w_2.rf");
      system.readWRGrid("out/testConversion1D_lam_w_2.rf");

      // Compare result to original
      RFieldComparison<1> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison3.maxDiff() < 1.0E-10);

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
      wFields_check = system.w().basis();

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
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      system.basisToKGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.kf");
      system.kGridToBasis("out/testConversion2D_hex_w.kf",
                          "out/testConversion2D_hex_w_2.bf");
      system.readWBasis("out/testConversion2D_hex_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.readWRGrid("out/testConversion2D_hex_w.rf");
      DArray< RField<2> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      system.rGridToKGrid("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w_2.kf");
      system.kGridToRGrid("out/testConversion2D_hex_w_2.kf",
                          "out/testConversion2D_hex_w_2.rf");
      system.readWRGrid("out/testConversion2D_hex_w_2.rf");

      // Compare result to original
      RFieldComparison<2> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff() << "\n";
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
      system.readWBasis("in/diblock/bcc/omega.in");

      // Store components of field as input
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

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
      BFieldComparison comparison1;
      comparison1.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison1.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison1.maxDiff() < 1.0E-10);

      // Round trip conversion basis -> kgrid -> basis, read result
      system.basisToKGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.kf");
      system.kGridToBasis("out/testConversion3D_bcc_w.kf",
                          "out/testConversion3D_bcc_w_2.bf");
      system.readWBasis("out/testConversion3D_bcc_w_2.bf");

      // Compare result to original
      BFieldComparison comparison2;
      comparison2.compare(wFields_check, system.w().basis());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison2.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison2.maxDiff() < 1.0E-10);

      // Round trip conversion rgrid -> kgrid -> rgrid, read result
      system.readWRGrid("out/testConversion3D_bcc_w.rf");
      DArray< RField<3> > wFieldsRGrid_check;
      wFieldsRGrid_check = system.w().rgrid();

      system.rGridToKGrid("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w_2.kf");
      system.kGridToRGrid("out/testConversion3D_bcc_w_2.kf",
                          "out/testConversion3D_bcc_w_2.rf");
      system.readWRGrid("out/testConversion3D_bcc_w_2.rf");

      // Compare result to original
      RFieldComparison<3> comparison3;
      comparison3.compare(wFieldsRGrid_check, system.w().rgrid());
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison3.maxDiff() << "\n";
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

      system.readWBasis("in/diblock/bcc/omega.in");
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
      wFields_check = system.w().basis();

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                     "out/testIterate1D_lam_rigid_w.bf", 
                                     system.w().basis(), 
                                     system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                     "out/testIterate1D_lam_rigid_c.bf", 
                                     system.c().basis(), 
                                     system.domain().unitCell());

      // Compare solution to original fields
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

      // Check stress value
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT(std::abs(stress) < 1.0E-8);

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

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                    "out/testIterate1D_lam_flex_w.bf", 
                                    system.w().basis(), 
                                    system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                    "out/testIterate1D_lam_flex_c.bf", 
                                    system.c().basis(), 
                                    system.domain().unitCell());

      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      system.readWBasis("in/diblock/lam/omega.ref");
      //wFields_check = system.w().basis();

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);
   }

   void testIterate1D_lam_soln()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_soln.log");

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/solution/lam/param", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/solution/lam/w.bf");
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate1D_lam_soln_w.bf", 
                                      system.w().basis(), 
                                      system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate1D_lam_soln_c.bf", 
                                      system.c().basis(), 
                                      system.domain().unitCell());
      system.writeBlockCRGrid("out/testIterate1D_lam_soln_block_c.rf");

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 2.0E-6);
   }

   void testIterate1D_lam_open_soln()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_soln.log");

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      
      std::ifstream in;
      openInputFile("in/solution/lam_open/param", in);
      system.readParam(in);
      in.close();

      // Read in comparison result
      system.readWBasis("in/solution/lam_open/w.ref");
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/solution/lam_open/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                   "out/testIterate1D_lam_open_soln_w.bf", 
                                   system.w().basis(), 
                                   system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                   "out/testIterate1D_lam_open_soln_c.bf", 
                                   system.c().basis(), 
                                   system.domain().unitCell());

      // Compare result
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
   }

   void testIterate1D_lam_open_blend()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_blend.log");

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      
      std::ifstream in;
      openInputFile("in/blend/lam/param", in);
      system.readParam(in);
      in.close();

      // Read in comparison result
      system.readWBasis("in/blend/lam/w.ref");
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/blend/lam/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                 "out/testIterate1D_lam_open_blend_w.bf", 
                                 system.w().basis(), 
                                 system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                 "out/testIterate1D_lam_open_blend_c.bf", 
                                 system.c().basis(), 
                                 system.domain().unitCell());
      system.writeBlockCRGrid("out/testIterate1D_lam_open_blend_block_c.rf");

      // Compare result
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-8);
   }

   void testIterate1D_lam_open_shift()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_shift.log");

      System<1> system, systemShift;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      systemShift.fileMaster().setInputPrefix(filePrefix());
      systemShift.fileMaster().setOutputPrefix(filePrefix());
      
      std::ifstream in, inShift;
      openInputFile("in/solution/lam_open/param", in);
      system.readParam(in);
      openInputFile("in/solution/lam_open/param", inShift);
      systemShift.readParam(inShift);
      in.close();

      // Read in comparison result
      system.readWBasis("in/solution/lam_open/w.ref");
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/solution/lam_open/w.bf");
      systemShift.readWBasis("in/solution/lam_open/w.bf");

      // Apply shift to input fields.
      double shift = 2;
      DArray<DArray <double> > wFields_ = systemShift.w().basis();
      for (int i = 0; i < systemShift.mixture().nMonomer(); ++i) {
         wFields_[i][0] += shift;
      }
      systemShift.setWBasis(wFields_);

      // Apply shift to polymer and solvent chemical potentials.
      for (int i = 0; i < systemShift.mixture().nSolvent(); ++i) {
         double L = systemShift.mixture().solvent(i).size();
         double newMu = systemShift.mixture().solvent(i).mu() + L*shift;
         systemShift.mixture().solvent(i).setMu(newMu);
      }
      for (int i = 0; i < systemShift.mixture().nPolymer(); ++i) {
         double L = 0;
         for (int j = 0; j < systemShift.mixture().polymer(i).nBlock(); ++j) {
            L += systemShift.mixture().polymer(i).block(j).length();
         }
         double newMu = systemShift.mixture().polymer(i).mu() + L*shift;
         systemShift.mixture().polymer(i).setMu(newMu);
      }

      // Iterate
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      int errorShift = systemShift.iterate();
      if (errorShift) {
         TEST_THROW("Shifted iterator failed to converge.");
      }

      // Verify concentration fields, thermo, and pressure
      BFieldComparison comparison(1);
      comparison.compare(system.c().basis(),systemShift.c().basis());
      double fDiff, pDiff;
      if (!system.hasFreeEnergy()) system.computeFreeEnergy();
      if (!systemShift.hasFreeEnergy()) systemShift.computeFreeEnergy();
      fDiff = std::abs(system.fHelmholtz() - systemShift.fHelmholtz());
      pDiff = std::abs(system.pressure() - systemShift.pressure() + shift);
      TEST_ASSERT(comparison.maxDiff() < 5.0E-8);
      TEST_ASSERT(fDiff < 1E-6);
      TEST_ASSERT(pDiff < 1E-6);

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
      wFields_check = system.w().basis();

      // Iterate, output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                       "out/testIterate2D_hex_rigid_w.bf", 
                                       system.w().basis(), 
                                       system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                       "out/testIterate2D_hex_rigid_c.bf", 
                                       system.c().basis(), 
                                       system.domain().unitCell());

      // Compare current solution to reference solution
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum error of 2.608E-7 occurs for the first star

      // Check stress
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT (std::abs(stress) < 1.0E-8);
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
      wFields_check = system.w().basis();

      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate2D_hex_flex_w.bf", 
                                      system.w().basis(), 
                                      system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate2D_hex_flex_c.bf", 
                                      system.c().basis(), 
                                      system.domain().unitCell());

      // Compare solution to reference solution
      BFieldComparison comparison(1);
      // setVerbose(1);
      comparison.compare(wFields_check, system.w().basis());
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
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

      // Read initial guess
      system.readWBasis("in/diblock/bcc/omega.ref");

      // Save copy of initial fields
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate3D_bcc_rigid_w.bf", 
                                      system.w().basis(), 
                                      system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                      "out/testIterate3D_bcc_rigid_c.bf", 
                                      system.c().basis(), 
                                      system.domain().unitCell());

      // Compare solution to reference solution
      BFieldComparison comparison(1); // Constructor argument 1 skips star 0
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum difference of 1.023E-7 occurs for the second star

      // Test that stress is small
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT(std::abs(stress) < 1.0E-7);

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
      wFields_check = system.w().basis();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                     "out/testIterate3D_bcc_flex_w.bf", 
                                     system.w().basis(), 
                                     system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis(
                                     "out/testIterate3D_bcc_flex_c.bf", 
                                     system.c().basis(), 
                                     system.domain().unitCell());

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum difference of 1.09288E-7 occurs for the second star

   }

   void testIterate3D_altGyr_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_altGyr_flex.log");

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      // Read parameter file
      std::ifstream in;
      openInputFile("in/triblock/altGyr/param", in);
      system.readParam(in);
      in.close();

      // Input a converged solution from PSCF Fortran
      system.readWBasis("in/triblock/altGyr/w.bf");

      // Make copy of input fields for later comparison
      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis("out/testIterate3D_altGyr_flex_w.bf", 
                                        system.w().basis(), 
                                        system.domain().unitCell());
      system.domain().fieldIo().writeFieldsBasis("out/testIterate3D_altGyr_flex_c.bf", 
                                        system.c().basis(), 
                                        system.domain().unitCell());
      system.writeBlockCRGrid("out/testIterate3D_altGyr_flex_block_c.rf");

      // Compare w fields
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-6);

      // Compare Helmoltz free energies
      if (!system.hasFreeEnergy()) system.computeFreeEnergy();
      double fHelmholtz = system.fHelmholtz();
      double fHelmholtzRef = 3.9642295402;     // from PSCF Fortran
      double fDiff = fHelmholtz - fHelmholtzRef;
      if (verbose() > 0) {
         std::cout << "fHelmholtz diff = " << fDiff << "\n";
      }
      TEST_ASSERT(std::abs(fDiff) < 1.0E-7);

      // Compare relaxed unit cell parameters
      double cellParam = system.domain().unitCell().parameter(0); 
      double cellParamRef = 2.2348701424;     // from PSCF Fortran
      double cellDiff = cellParam - cellParamRef;
      if (verbose() > 0) {
         std::cout << "Cell param diff = " << cellDiff << "\n";
      }
      TEST_ASSERT(std::abs(cellDiff) < 1.0E-7);
   }

   void testIterate3D_c15_1_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_c15_1_flex.log");

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/diblock/c15_1/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/c15_1/w_ref.bf");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.w().basis();

      system.readWBasis("in/diblock/c15_1/w_in.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.domain().fieldIo().writeFieldsBasis(
                                   "out/testIterate3D_c15_1_flex_w.bf", 
                                    system.w().basis(), 
                                    system.domain().unitCell());
      system.domain().fieldIo().writeFieldsRGrid(
                                    "out/testIterate3D_c15_1_flex_w.rf", 
                                    system.c().rgrid(), 
                                    system.domain().unitCell());

      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.w().basis());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum difference of 1.09288E-7 occurs for the second star

   }

   void testIterateWithMaskAndH() // test manual entry of mask and h fields
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterateWithMaskAndH.log");
      
      // Set up system
      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/maskAndH/param", in);
      system.readParam(in);
      in.close();

      // Read initial guess
      system.readWBasis("in/maskAndH/w.bf");

      // Read in the mask and external fields from file
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      unitCell = system.domain().unitCell();
      system.mask().setFieldIo(system.domain().fieldIo());
      system.mask().allocate(system.domain().basis().nBasis(), 
                             system.domain().mesh().dimensions());
      system.mask().readBasis("in/maskAndH/mask.bf", unitCell);
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      system.h().setFieldIo(system.domain().fieldIo());
      system.h().allocateBasis(system.domain().basis().nBasis());
      system.h().allocateRGrid(system.domain().mesh().dimensions());
      system.h().readBasis("in/maskAndH/h.bf", unitCell);

      // Run the solve function
      system.iterate();

      // Check converged field is correct by comparing to files in in/maskAndH
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/maskAndH/w.ref", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);

      double diff = bComparison.maxDiff();
      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
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
TEST_ADD(SystemTest, testIterate1D_lam_soln)
TEST_ADD(SystemTest, testIterate1D_lam_open_soln)
TEST_ADD(SystemTest, testIterate1D_lam_open_blend)
TEST_ADD(SystemTest, testIterate1D_lam_open_shift)
TEST_ADD(SystemTest, testIterate2D_hex_rigid)
TEST_ADD(SystemTest, testIterate2D_hex_flex)
TEST_ADD(SystemTest, testIterate3D_bcc_rigid)
TEST_ADD(SystemTest, testIterate3D_bcc_flex)
TEST_ADD(SystemTest, testIterate3D_altGyr_flex)
TEST_ADD(SystemTest, testIterate3D_c15_1_flex)
TEST_ADD(SystemTest, testIterateWithMaskAndH)
TEST_END(SystemTest)

#endif
