#ifndef RPC_AM_ITERATOR_BASIS_TEST_H
#define RPC_AM_ITERATOR_BASIS_TEST_H

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
   * Initialize system - read parameter file and w-field basis file. 
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
      outFileRoot = makeFileRoot("out/testIterateBasis", outSuffix, D);

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
      //setVerbose(1);

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

      // Unit cell parameter
      double a  = system.domain().unitCell().parameter(0);
      if (verbose() > 0) {
         std::cout << "\na = " << Dbl(a, 20, 12);
      }
      TEST_ASSERT(std::abs(a - 1.5161093077) < 1.0E-8);

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
      system.computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 3.3354925e-01)/s0) < 1.0E-5);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.scft().fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.scft().fHelmholtz();

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
      system.computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 3.3354725E-3)/s0) < 1.0E-4);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.scft().fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.scft().fHelmholtz();

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
      if (verbose() > 0) {
         std::cout << "\n wMaxDiff = " << wMaxDiff;
         std::cout << "\n cMaxDiff = " << cMaxDiff;
      }
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
                  "in/blend/lam/param.open",
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

      // Process system
      System<1> system;
      initSystem(system,
                 "in/solution/lam_open/param",
                 "in/solution/lam_open/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }

      // Initialize systemShift
      System<1> systemShift;
      initSystem(systemShift,
                 "in/solution/lam_open/param",
                 "in/solution/lam_open/w.bf");

      // Shift w fields in SystemShift
      double shift = 2.0;
      DArray<DArray <double> > wFields_ = systemShift.w().basis();
      for (int i = 0; i < systemShift.mixture().nMonomer(); ++i) {
         wFields_[i][0] += shift;
      }
      systemShift.w().setBasis(wFields_);

      // Shift chemical potentials in SystemShift
      double L, newMu;
      int nSolvent = systemShift.mixture().nSolvent();
      for (int i = 0; i < nSolvent; ++i) {
         L = systemShift.mixture().solvent(i).size();
         newMu = systemShift.mixture().solvent(i).mu() + L*shift;
         systemShift.mixtureModifier().setMuSolvent(i, newMu);
      }
      int nPolymer = systemShift.mixture().nPolymer();
      int nBlock;
      for (int i = 0; i < nPolymer; ++i) {
         L = 0.0;
         nBlock = systemShift.mixture().polymer(i).nBlock();
         for (int j = 0; j < nBlock; ++j) {
            L += systemShift.mixture().polymer(i).block(j).length();
         }
         newMu = systemShift.mixture().polymer(i).mu() + L*shift;
         systemShift.mixtureModifier().setMuPolymer(i, newMu);
      }

      // Iterate systemShift
      int errorShift = systemShift.iterate();
      if (errorShift) {
         TEST_THROW("Shifted iterator failed to converge.");
      }

      // Compare converged concentration fields
      BFieldComparison comparison(1);
      comparison.compare(system.c().basis(), systemShift.c().basis());
      TEST_ASSERT(comparison.maxDiff() < 5.0E-8);

      // Compare free energies and pressures
      double fDiff, pDiff;
      if (!system.scft().hasData()) system.scft().compute();
      if (!systemShift.scft().hasData()) systemShift.scft().compute();
      fDiff = std::abs(system.scft().fHelmholtz() - systemShift.scft().fHelmholtz());
      pDiff = std::abs(system.scft().pressure() - systemShift.scft().pressure() + shift);
      TEST_ASSERT(fDiff < 1.0E-6);
      TEST_ASSERT(pDiff < 1.0E-6);
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

      FieldIo<2> const & fieldIo = system.domain().fieldIo();
      fieldIo.scaleFieldsBasis("out/testIterateBasis2D_hex_flex_w.bf",
                               "out/testIterateBasis2D_hex_flex_w_scaled.bf",
                               0.01);
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
      system.computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 1.6546466e-01)/s0) < 1.0E-5);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.scft().fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.scft().fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 1.0E-4);
      #if 0
      #endif

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
      outFileRoot = makeFileRoot("out/testIterateBasis", "hex_bead_stress", 2);
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
      system.computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs((s0 - 1.654755884391e-03)/s0) < 1.0E-5);

      // Solve for lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.scft().fHelmholtz();

      // Solve for upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.scft().fHelmholtz();

      // Numerical stress by finite differences
      double sn = (fp - fm)/(ap - am);
      if (verbose() > 0) {
         std::cout << "\nStress [scft]  = " << Dbl(s0, 20, 12);
         std::cout << "\nStress [diff]  = " << Dbl(sn, 20, 12);
         std::cout << "\nError fraction = " << Dbl((s0-sn)/s0, 20, 12);
      }
      TEST_ASSERT(std::abs((sn - s0)/s0) < 4.0E-4);

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

      // Compare unit cell value
      double a = system.domain().unitCell().parameter(0);
      //std::cout << "a = " << Dbl(a,20,12) << std::endl;
      TEST_ASSERT(std::abs(a - 1.7593559883) < 1.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
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
      system.computeStress();
      double s0 = system.mixture().stress(0);
      TEST_ASSERT(std::fabs(s0 - 8.0561073E-2) < 1.0E-6);

      // Lower unit cell parameter am (minus)
      parameters.clear();
      parameters.append(am);
      system.setUnitCell(parameters);
      system.iterate();
      double fm = system.scft().fHelmholtz();
      
      // Upper unit cell parameter ap (plus)
      parameters.clear();
      parameters.append(ap);
      system.setUnitCell(parameters);
      system.iterate();
      double fp = system.scft().fHelmholtz();

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

      system.mixture().writeBlockCRGrid("out/testIterate3D_altGyr_flex_block_c.rf");

      // Compare Helmoltz free energies
      if (!system.scft().hasData()) system.scft().compute();
      double fHelmholtz = system.scft().fHelmholtz();
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

      // Compare to input reference w_ref.bf
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
      initSystem(system,"in/maskAndH/param", "in/maskAndH/w.bf");

      // Read in the mask and external fields from file
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      unitCell = system.domain().unitCell();
      system.mask().setFieldIo(system.domain().fieldIo());
      system.mask().allocateBasis(system.domain().basis().nBasis());
      system.mask().allocateRGrid(system.domain().mesh().dimensions());
      system.mask().readBasis("in/maskAndH/mask.bf");
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      system.h().setFieldIo(system.domain().fieldIo());
      system.h().allocateBasis(system.domain().basis().nBasis());
      system.h().allocateRGrid(system.domain().mesh().dimensions());
      system.h().readBasis("in/maskAndH/h.bf");

      // Iterate to convergence
      system.iterate();

      // Compare converged field to in/maskAndH/w.ref
      double diff = readCompareWBasis(system, "in/maskAndH/w.ref");
      TEST_ASSERT(diff < 1.0E-5);
   }

};

TEST_BEGIN(AmIteratorBasisTest)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_flex)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_stress)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_bead_flex)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_bead_stress)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_soln)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_soln)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_blend)
TEST_ADD(AmIteratorBasisTest, testIterate1D_lam_open_shift)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_flex)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_stress)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_bead_flex)
TEST_ADD(AmIteratorBasisTest, testIterate2D_hex_bead_stress)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_rigid)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_bcc_stress)
TEST_ADD(AmIteratorBasisTest, testIterate3D_altGyr_flex)
TEST_ADD(AmIteratorBasisTest, testIterate3D_c15_1_flex)
TEST_ADD(AmIteratorBasisTest, testIterateWithMaskAndH)
TEST_END(AmIteratorBasisTest)

#endif
