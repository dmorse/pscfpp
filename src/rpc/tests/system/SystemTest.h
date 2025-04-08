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

   void tearDown()
   {  setVerbose(0); }

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;
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
   * Read basis fields into an array.
   */
   double compareBasis(DArray< DArray<double> > const & fields1,
                       DArray< DArray<double> > const & fields2)
   {
      BFieldComparison comparison(1);
      comparison.compare(fields1, fields2);
      double maxDiff = comparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "max error = " << maxDiff;
      }
      return maxDiff;
   }

   /*
   * Read basis fields into an array.
   */
   template <int D>
   double readCompareWBasis(System<D> const & system,
                            std::string filename)
   {
      DArray< DArray<double> > fields;
      UnitCell<D> unitCell;
      readBasisFields(system, filename, fields, unitCell);
      double maxDiff = compareBasis(fields, system.w().basis());
      return maxDiff;
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
      UTIL_CHECK(system.hasFreeEnergy());
      fHelmholtz = std::abs(fHelmholtz - system.fHelmholtz());
      pressure   = std::abs(pressure - system.pressure());
      if (verbose() > 0) {
         std::cout << "\nfHelmholtz diff = " << fHelmholtz;
         std::cout << "\npressure diff   = " << pressure;
      }
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
                    bool compareCFields = true,
                    bool compareToInput = false)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::string outFileRoot = "out/testIterate";
      outFileRoot += std::to_string(D) + "D_" + outSuffix;

      // Open log file
      std::string filename;
      filename = outFileRoot + ".log";
      openLogFile(filename.c_str());

      // Read parameter file
      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis(wFileName);

      // Construct reference fields for comparison
      DArray< DArray<double> > w_ref;
      DArray< DArray<double> > c_ref;
      UnitCell<D> unitCell;
      if (compareToInput) {
         w_ref = system.w().basis();
      } else {
         filename = "ref/";
         filename += "testIterate" + std::to_string(D) + "D_" + outSuffix;
         readBasisFields(system, filename + "_w.bf", w_ref, unitCell);
         readBasisFields(system, filename + "_c.bf", c_ref, unitCell);
      }

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      FieldIo<D> const & fieldIo = system.domain().fieldIo();
      fieldIo.writeFieldsBasis(outFileRoot + "_w.bf", 
                               system.w().basis(), 
                               system.domain().unitCell());
      fieldIo.writeFieldsBasis(outFileRoot + "_c.bf", 
                               system.c().basis(), 
                               system.domain().unitCell());

      // Compare solution to reference w fields
      BFieldComparison comparison(1);
      comparison.compare(w_ref, system.w().basis());
      wMaxDiff = comparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "w field max error = " << wMaxDiff;
      }

      // Compare solution to reference c fields
      if (compareCFields) {
         BFieldComparison comparison(1);
         comparison.compare(c_ref, system.c().basis());
         cMaxDiff = comparison.maxDiff();
         if (verbose() > 0) {
            std::cout << "\n";
            std::cout << "c field max error = " << cMaxDiff;
         }
      }

   }

   template <int D>
   FSArray<double, 6> computeStress(System<D>& system)
   {
      FSArray<double, 6> stress;
      int nParameter = system.domain().unitCell().nParameter();

      system.mixture().computeStress();
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
         std::cout << "Max error = " << comparison1.maxDiff();
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
         std::cout << "Max error = " << comparison2.maxDiff();
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
         std::cout << "Max error = " << comparison3.maxDiff();
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
         std::cout << "Max error = " << comparison1.maxDiff();
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
         std::cout << "Max error = " << comparison2.maxDiff();
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
         std::cout << "Max error = " << comparison1.maxDiff();
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
         std::cout << "Max error = " << comparison2.maxDiff();
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

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/lam/param.rigid",
                  "in/diblock/lam/omega.ref",
                  "lam_rigid",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare to input field omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/lam/omega.ref");
      TEST_ASSERT(wMaxDiff < 2.0E-7);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.42932542391e+00;  // Value from v1.1
      double pressure = 3.01212693885e+00;    // Value from v1.1
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.ref as input, and compared to input

   }

   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/lam/param.flex",
                  "in/diblock/lam/omega.in",
                  "lam_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 5.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =  2.42932542391e+00;
      double pressure   =  3.01212693880e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to reference omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/lam/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   void testIterate1D_lam_soln()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/solution/lam/param",
                  "in/solution/lam/w.bf",
                  "lam_soln",
                  wMaxDiff,
                  cMaxDiff);
      //std::cout << "\n wMaxDiff = " << wMaxDiff;
      //std::cout << "\n cMaxDiff = " << cMaxDiff;
      TEST_ASSERT(wMaxDiff < 5.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =  -2.75154924314e+01;
      double pressure   =  3.24415250702e+01;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input w.bf
      wMaxDiff = readCompareWBasis(system, "in/solution/lam/w.bf");
      TEST_ASSERT(wMaxDiff < 2.0E-6);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to input
   }

   void testIterate1D_lam_open_soln()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/solution/lam_open/param",
                  "in/solution/lam_open/w.bf",
                  "lam_open_soln",
                  wMaxDiff,
                  cMaxDiff);
      //std::cout << "\n wMaxDiff = " << wMaxDiff;
      //std::cout << "\n cMaxDiff = " << cMaxDiff;
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = -2.75154924224e+01;
      double pressure   =  3.24415250699e+01;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to reference w.ref
      wMaxDiff = readCompareWBasis(system, "in/solution/lam_open/w.ref");
      TEST_ASSERT(wMaxDiff < 2.0E-6);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to w.ref

   }

   void testIterate1D_lam_open_blend()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/blend/lam/param",
                  "in/blend/lam/w.bf",
                  "lam_open_blend",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = -3.30485085194e+00;
      double pressure = 5.10158598280e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to reference w.ref
      wMaxDiff = readCompareWBasis(system, "in/blend/lam/w.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to w.ref
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

      System<2> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/hex/param.rigid",
                  "in/diblock/hex/omega.ref",
                  "hex_rigid",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.80222103795e+00;
      double pressure = 3.19716573940e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/hex/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);
      // Maximum error of 2.608E-7 occurs for the first star

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.ref as input, compared to input
   }

   void testIterate2D_hex_flex()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/hex/param.flex",
                  "in/diblock/hex/omega.in",
                  "hex_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   2.80222103517e+00;
      double pressure =     3.19716584873e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/hex/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);
      // Maximum difference of 2.58E-7 occurs for the first star

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/bcc/param.rigid",
                  "in/diblock/bcc/omega.ref",
                  "bcc_rigid",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918380842e+00;
      double pressure =     4.03176984344e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/bcc/omega.ref");
      TEST_ASSERT(wMaxDiff < 1.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.ref as input, compared to input
   }

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/bcc/param.flex",
                  "in/diblock/bcc/omega.in",
                  "bcc_flex",
                  wMaxDiff,
                  cMaxDiff);
      //std::cout << "\n wMaxDiff = " << wMaxDiff;
      //std::cout << "\n cMaxDiff = " << cMaxDiff;
      TEST_ASSERT(wMaxDiff < 1.0E-6);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918376624e+00;
      double pressure =     4.03176988267e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/bcc/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   void testIterate3D_altGyr_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/triblock/altGyr/param",
                  "in/triblock/altGyr/w.bf",
                  "altGyr_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare to input w.bf
      wMaxDiff = readCompareWBasis(system, "in/triblock/altGyr/w.bf");
      TEST_ASSERT(wMaxDiff < 5.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to input

      system.writeBlockCRGrid("out/testIterate3D_altGyr_flex_block_c.rf");

      // Compare Helmoltz free energies
      if (!system.hasFreeEnergy()) system.computeFreeEnergy();
      double fHelmholtz = system.fHelmholtz();
      double fHelmholtzRef = 3.9642295402;     // from PSCF Fortran
      double fDiff = fHelmholtz - fHelmholtzRef;
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "fHelmholtz diff = " << fDiff;
      }
      TEST_ASSERT(std::abs(fDiff) < 1.0E-7);

      // Compare relaxed unit cell parameters
      double cellParam = system.domain().unitCell().parameter(0); 
      double cellParamRef = 2.2348701424;     // from PSCF Fortran
      double cellDiff = cellParam - cellParamRef;
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell param diff = " << cellDiff;
      }
      TEST_ASSERT(std::abs(cellDiff) < 1.0E-7);
   }

   void testIterate3D_c15_1_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/c15_1/param.flex",
                  "in/diblock/c15_1/w_in.bf",
                  "c15_1_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-6);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.21902602858e+00;
      double pressure =     3.37292021160e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input w.bf
      wMaxDiff = readCompareWBasis(system, "in/diblock/c15_1/w_ref.bf");
      TEST_ASSERT(wMaxDiff < 1.0E-6);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-7);

      // v1.1 test used w_in.bf as input, compares to w_ref.bf
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
      system.mask().allocateBasis(system.domain().basis().nBasis()); 
      system.mask().allocateRGrid(system.domain().mesh().dimensions());
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
