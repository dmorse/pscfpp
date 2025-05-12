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
   * Read basis fields into an array.
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
   * Read basis fields into an array.
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
   * Template for an iteration test, with regression testing.
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
      system.readWBasis(wFileName);
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
      system.writeWBasis(outFileRoot + "_w.bf");
      system.writeCBasis(outFileRoot + "_w.bf");
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
      outFileRoot = makeFileRoot("out/testIterateBasis", outSuffix, D);

      // Open log file
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

   void testIterate1D_lam_stress()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateBasis", "lam_stress", 1);
      openLogFile(outFileRoot + ".log");

      // Initialize, read reference solution
      System<1> system;
      initSystem(system, 
                 "in/diblock/lam/param.rigid",
                 "ref/testIterate1D_lam_flex_w.bf");
      system.iterate();
      double ar  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(ar,20, 12);
      }

      // Unit cell size values
      double a0  = ar * 1.200;
      double am  = a0 * 0.995;
      double ap  = a0 * 1.005;

      FSArray<double, 6> parameters;
      parameters.clear();

      // Solve for middle unit cell parameter a0 
      parameters.append(a0);
      system.setUnitCell(parameters);
      system.iterate();
      system.mixture().computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 3.3354925e-01)/s0) < 1.0E-5);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 1.0E-4);

   }

   void testIterate1D_lam_bead_flex()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      testIterate(system,
                  "in/diblock/lam_bead/param.flex",
                  "in/diblock/lam_bead/w.bf",
                  "lam_bead_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 5.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Unit cell parameter
      double a  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(a, 20, 12);
      }
      TEST_ASSERT(std::abs(a - 1.51567153218) < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.42834917413e-02;
      double pressure   = 3.00896311324e-02;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-8);
      TEST_ASSERT(pressure < 1.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-10);
   }

   void testIterate1D_lam_bead_stress()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateBasis", "lam_bead_stress", 1);
      openLogFile(outFileRoot + ".log");

      // Initialize, read reference solution
      System<1> system;
      initSystem(system, 
                 "in/diblock/lam_bead/param.rigid",
                 "ref/testIterate1D_lam_bead_flex_w.bf");
      system.iterate(); 
      double ar  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(ar,20, 12);
      }

      // Unit cell size values
      double a0  = ar * 1.200;
      double am  = a0 * 0.995;
      double ap  = a0 * 1.005;

      FSArray<double, 6> parameters;
      parameters.clear();

      // Solve for middle unit cell parameter a0 
      parameters.append(a0);
      system.setUnitCell(parameters);
      system.iterate();
      system.mixture().computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 3.3354725E-3)/s0) < 1.0E-4);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 1.0E-4);

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

   void testIterate1D_lam_blend()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);
        
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


   void testIterate2D_hex_stress()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateBasis", "hex_stress", 2);
      openLogFile(outFileRoot + ".log");

      // Initialize, read reference solution
      System<2> system;
      initSystem(system, 
                 "in/diblock/hex/param.rigid",
                 "ref/testIterate2D_hex_flex_w.bf");
      system.iterate();
      double ar  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(ar, 20, 12);
      }

      // Unit cell size values
      double a0  = ar * 1.100;
      double am  = a0 * 0.996;
      double ap  = a0 * 1.004;

      FSArray<double, 6> parameters;
      parameters.clear();

      // Solve for middle unit cell parameter a0 
      parameters.append(a0);
      system.setUnitCell(parameters);
      system.iterate();
      system.mixture().computeStress();
      double s0 = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "\nstress = " << Dbl(s0, 20, 12);
      }
      TEST_ASSERT(std::fabs((s0 - 1.6546466e-01)/s0) < 1.0E-5);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 1.0E-4);

   }

   void testIterate2D_hex_bead_flex()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      double wMaxDiff, cMaxDiff;
      System<2> system;
      testIterate(system,
                  "in/diblock/hex_bead/param.flex",
                  "in/diblock/hex_bead/w.bf",
                  "hex_bead_flex",
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.80048919599e-02;
      double pressure =  3.19132478045e-02;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);
   }

   void testIterate2D_hex_bead_stress()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateBasis", 
                                 "hex_bead_stress", 2);
      openLogFile(outFileRoot + ".log");

      // Initialize, read reference solution
      System<2> system;
      initSystem(system, 
                 "in/diblock/hex_bead/param.rigid",
                 "ref/testIterate2D_hex_bead_flex_w.bf");
      system.iterate(); 
      double ar  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(ar,20, 12);
      }

      // Unit cell size values
      double a0  = ar * 1.100;
      double am  = a0 * 0.995;
      double ap  = a0 * 1.005;

      FSArray<double, 6> parameters;
      parameters.clear();

      // Solve for middle unit cell parameter a0 
      parameters.append(a0);
      system.setUnitCell(parameters);
      system.iterate();
      system.mixture().computeStress();
      double s0 = system.mixture().stress(0);
      double diff = std::fabs((s0 - 1.654755884391e-03)/s0);
      if (verbose() > 0) {
         std::cout << "\nstress = " << Dbl(s0, 20, 12);
         std::cout << "\ndiff   = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-6);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 4.0E-4);

   }


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

   void testIterate3D_bcc_stress()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateBasis", "bcc_stress", 3);
      openLogFile(outFileRoot + ".log");

      System<3> system;
      initSystem(system,
                 "in/diblock/bcc/param.rigid",
                 "in/diblock/bcc/omega.in");

      // Unit cell size values
      double ar  = system.domain().unitCell().parameter(0);
      double a0  = ar*1.100;
      double am  = a0*0.995;
      double ap  = a0*1.005;

      FSArray<double, 6> parameters;

      // Middle unit cell parameter a0 
      parameters.append(a0);
      system.setUnitCell(parameters);
      system.iterate();
      system.mixture().computeStress();
      double s0 = system.mixture().stress(0);
      double diff = std::fabs((s0 - 8.0561073E-2)/s0);
      if (verbose() > 0) {
	  std::cout << "\nstress  = " << s0;
	  std::cout << "\ndiff    = " << s0;
      }
      TEST_ASSERT(diff < 1.0E-5);

      // Lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.fHelmholtz();
      
      // Upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 1.0E-3);

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

TEST_BEGIN(AmIteratorBasisTest)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_flex)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_stress)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_bead_flex)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_bead_stress)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_soln)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_blend)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_blend)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_soln)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_flex)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_stress)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_bead_flex)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_bead_stress)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_flex_noBatched)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_stress)
TEST_ADD(AmIteratorBasisTest, testIterate3D_altGyr_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_c15_1_flex)
TEST_ADD(AmIteratorBasisTest, testIterateWithMaskAndH)

TEST_END(AmIteratorBasisTest)

#endif
