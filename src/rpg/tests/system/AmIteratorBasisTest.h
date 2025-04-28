#ifndef RPG_AM_ITERATOR_TEST_H
#define RPG_AM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>
#include <prdc/crystal/BFieldComparison.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class AmIteratorBasisTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

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

   template <int D>
   void setupSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
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

      std::string outFileRoot = "out/testIterateBasis";
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
      system.writeWBasis(outFileRoot + "_w.bf");
      system.writeCBasis(outFileRoot + "_c.bf");
      //FieldIo<D> const & fieldIo = system.domain().fieldIo();
      //fieldIo.writeFieldsBasis(outFileRoot + "_w.bf",
      //                         system.w().basis(),
      //                         system.domain().unitCell());
      //fieldIo.writeFieldsBasis(outFileRoot + "_c.bf",
      //                         system.c().basis(),
      //                         system.domain().unitCell());

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

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.42932542391e+00;  // Value from v1.1
      double pressure = 3.01212693885e+00;    // Value from v1.1
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input field omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/lam/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.ref as input, and compared to input
   }

   #if 0
   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_rigid.log");

      System<1> system;
      setupSystem<1>(system,"in/diblock/lam/param.rigid");

      // Read w fields
      system.readWBasis("in/diblock/lam/omega.ref");

      // Make reference copy of w fields in basis format
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate1D_lam_rigid_w.bf");
      system.writeCBasis("out/testIterate1D_lam_rigid_c.bf");

      // Compare result to original in basis format
      BFieldComparison comparison(1);
      comparison.compare(b_wFields_check, system.w().basis());
      if (verbose() > 0) {
         Log::file() << "\n";
         Log::file() << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
   }
   #endif

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

   #if 0
   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_flex.log");

      System<1> system;
      setupSystem<1>(system,"in/diblock/lam/param.flex");

      system.readWBasis("in/diblock/lam/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_flex_w.bf");
      system.writeCBasis("out/testIterate1D_lam_flex_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (verbose() > 0 || diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate1D_lam_flex_noBatched()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterateBasis1D_lam_flex_noBatched.log");

      System<1> system;
      setupSystem<1>(system,"in/diblock/lam/param_noBatched.flex");

      system.readWBasis("in/diblock/lam/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterateBasis1D_lam_flex_noBatched_w.bf");
      system.writeCBasis("out/testIterateBasis1D_lam_flex_noBatched_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (verbose() > 0 || diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
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

   #if 0
   void testIterate1D_lam_soln()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_soln.log");

      System<1> system;
      setupSystem<1>(system,"in/solution/lam/param");

      // Make reference copy of w fields
      system.readWBasis("in/solution/lam/w.bf");
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_soln_w.bf");
      system.writeCBasis("out/testIterate1D_lam_soln_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 2.0E-6;
      if (verbose() > 0 || diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate1D_lam_blend()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterateBasis1D_lam_blend.log");

      System<1> system;
      setupSystem<1>(system,"in/blend/lam/param.closed");

      // Make reference copy of w fields
      system.readWBasis("in/blend/lam/w.ref");
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/blend/lam/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterateBasis1D_lam_blend_w.bf");
      system.writeCBasis("out/testIterateBasis1D_lam_blend_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      //setVerbose(1);
      if (verbose() > 0 || diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
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
      TEST_ASSERT(wMaxDiff < 5.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare to reference w.ref
      wMaxDiff = readCompareWBasis(system, "in/solution/lam_open/w.ref");
      TEST_ASSERT(wMaxDiff < 2.0E-6);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to w.ref

   }

   #if 0
   void testIterate1D_lam_open_soln()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_soln.log");

      System<1> system;
      setupSystem<1>(system,"in/solution/lam_open/param");

      // Make reference copy of w fields
      system.readWBasis("in/solution/lam_open/w.ref");
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/solution/lam_open/w.bf");

      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_open_soln_w.bf");
      system.writeCBasis("out/testIterate1D_lam_open_soln_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      //double epsilon = 5.0E-7;
      double epsilon = 6.0E-6;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate1D_lam_open_blend()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/blend/lam/param.open",
                  "in/blend/lam/w.bf",
                  "lam_open_blend",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare to reference w.ref
      wMaxDiff = readCompareWBasis(system, "in/blend/lam/w.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to w.ref
   }

   #if 0
   void testIterate1D_lam_open_blend()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_blend.log");

      System<1> system;
      setupSystem<1>(system,"in/blend/lam/param.open");

      // Make reference copy of w fields
      system.readWBasis("in/blend/lam/w.ref");
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/blend/lam/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_open_blend_w.bf");
      system.writeCBasis("out/testIterate1D_lam_open_blend_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/hex/param.rigid",
                  "in/diblock/hex/omega.ref",
                  "hex_rigid", wMaxDiff, cMaxDiff);
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

   #if 0
   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_rigid.log");

      System<2> system;
      setupSystem<2>(system,"in/diblock/hex/param.rigid");

      // Read reference solution
      system.readWBasis("in/diblock/hex/omega.ref");
      TEST_ASSERT(system.domain().basis().isInitialized());

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      // Read initial guess, iterate, output solution
      // RPC tests start from the reference solution,
      // rather than a nearby solution, so I guess do that here too?
      // system.readWBasis("in/diblock/hex/omega.in");
      Log::file() << "Beginning iteration \n";
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate2D_hex_rigid_w.bf");
      system.writeCBasis("out/testIterate2D_hex_rigid_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate2D_hex_flex()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/hex/param.flex",
                  "in/diblock/hex/omega.in",
                  "hex_flex", wMaxDiff, cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-6);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   2.80222103517e+00;
      double pressure =     3.19716584873e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/hex/omega.ref");
      TEST_ASSERT(wMaxDiff < 1.0E-6);
      // Maximum difference of 2.58E-7 occurs for the first star

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   #if 0
   void testIterate2D_hex_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_flex.log");

      System<2> system;
      setupSystem<2>(system,"in/diblock/hex/param.flex");

      // Read reference solution (produced by Fortran code)
      system.readWBasis("in/diblock/hex/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate2D_hex_flex_w.bf");
      system.writeCBasis("out/testIterate2D_hex_flex_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 8.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate2D_hex_flex_noBatched()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterateBasis2D_hex_flex_noBatched.log");

      System<2> system;
      setupSystem<2>(system,"in/diblock/hex/param_noBatched.flex");

      // Read reference solution (produced by Fortran code)
      system.readWBasis("in/diblock/hex/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterateBasis2D_hex_flex_noBatched_w.bf");
      system.writeCBasis("out/testIterateBasis2D_hex_flex_noBatched_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 8.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/bcc/param.rigid",
                  "in/diblock/bcc/omega.ref",
                  "bcc_rigid", wMaxDiff, cMaxDiff);
      TEST_ASSERT(wMaxDiff < 5.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-7);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918380842e+00;
      double pressure =     4.03176984344e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to input omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/bcc/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.ref as input, compared to input
   }

   #if 0
   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_rigid.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param.rigid");

      system.readWBasis("in/diblock/bcc/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_rigid_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_rigid_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 7.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/bcc/param.flex", 
		  "in/diblock/bcc/omega.in",
                  "bcc_flex", wMaxDiff, cMaxDiff);
      TEST_ASSERT(wMaxDiff < 5.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918376624e+00;
      double pressure =     4.03176988267e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare to reference omega.ref
      wMaxDiff = readCompareWBasis(system, "in/diblock/bcc/omega.ref");
      TEST_ASSERT(wMaxDiff < 5.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   #if 0
   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_flex.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param.flex");

      system.readWBasis("in/diblock/bcc/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_flex_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_flex_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
   }
   #endif

   void testIterate3D_bcc_flex_noBatched()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterateBasis3D_bcc_flex_noBatched.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param_noBatched.flex");

      system.readWBasis("in/diblock/bcc/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      b_wFields_check = system.w().basis();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterateBasis3D_bcc_flex_noBatched_w.bf");
      system.writeCBasis("out/testIterateBasis3D_bcc_flex_noBatched_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);

      // Compare difference to tolerance epsilon
      double diff = comparison.maxDiff();
      double epsilon = 5.0E-7;
      if (diff > epsilon) {
         Log::file() << "\n";
         Log::file() << "Max diff = " << comparison.maxDiff() << "\n";
         Log::file() << "Rms diff = " << comparison.rmsDiff() << "\n";
         Log::file() << "epsilon  = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);
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

      system.writeBlockCRGrid("out/testIterateBasis3D_altGyr_flex_block_c.rf");

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
      unitCell = system.unitCell();
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

TEST_BEGIN(AmIteratorBasisTest)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_flex)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_soln)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_blend)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_blend)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_soln)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_flex)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate3D_altGyr_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_c15_1_flex)
TEST_ADD(AmIteratorBasisTest, testIterateWithMaskAndH)

TEST_END(AmIteratorBasisTest)

#endif
