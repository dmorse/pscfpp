#ifndef RPC_AM_ITERATOR_GRID_TEST_H
#define RPC_AM_ITERATOR_GRID_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/system/System.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldComparison.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpc;

class AmIteratorGridTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {  
      setVerbose(0); 
      closeLogFile();
   }

   // Utility functions

   /*
   * Allocate an array of rgrid fields.
   */
   template <int D>
   void allocateRGridFields(System<D> const & system,
                            DArray< RField<D> >& fields)
   {
      // Check and allocate outer DArray
      int nMonomer = system.mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      if (!fields.isAllocated()) {
         fields.allocate(nMonomer);
      }
      UTIL_CHECK(fields.capacity() == nMonomer);

      // Allocate fields
      Mesh<D> const & mesh = system.domain().mesh();
      IntVec<D> const & meshDimensions = mesh.dimensions();
      int meshSize = mesh.size();
      UTIL_CHECK(meshSize > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (!fields[i].isAllocated()) {
            fields[i].allocate(meshDimensions);
         }
         UTIL_CHECK(fields[i].capacity() == meshSize);
      }

   }

   /*
   * Read r-grid fields into an array.
   */
   template <int D>
   void readRGridFields(System<D> const & system,
                        std::string filename,
                        DArray< RField<D> >& fields,
                        UnitCell<D>& unitCell)
   {
      allocateRGridFields(system, fields);
      FieldIo<D> const & fieldIo = system.domain().fieldIo();
      fieldIo.readFieldsRGrid(filename, fields, unitCell);
   }

   /*
   * Compare rgrid fields.
   */
   template <int D>
   double compareRGrid(DArray< RField<D> > const & fields1,
                       DArray< RField<D> > const & fields2,
                       std::string message)
   {
      RFieldComparison<D> comparison;
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
   double readCompareWRGrid(System<D> const & system,
                            std::string filename)
   {
      DArray< RField<D> > fields;
      UnitCell<D> unitCell;
      readRGridFields(system, filename, fields, unitCell);
      std::string message = "max w field error = ";
      return compareRGrid(fields, system.w().rgrid(), message);
   }

   /*
   * Compare system w fields to reference fields from a file.
   */
   template <int D>
   double readCompareCRGrid(System<D> const & system,
                            std::string filename)
   {
      DArray< RField<D> > fields;
      UnitCell<D> unitCell;
      readRGridFields(system, filename, fields, unitCell);
      std::string message = "max c field error = ";
      return compareRGrid(fields, system.c().rgrid(), message);
   }

   template <int D>
   void setupSystem(System<D>& system, std::string paramFileName)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(paramFileName, in);
      system.readParam(in);
      in.close();
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
   * Template for an iteration test, with regression testing.
   */
   template <int D>
   void testIterate(System<D>& system,
                    std::string paramFileName,
                    std::string wFileName,
                    std::string outSuffix,
                    int& error,
                    double& wMaxDiff,
                    double& cMaxDiff,
                    bool compareCFields = true,
                    bool convertRef = true)
   {
      std::string outFileRoot;
      outFileRoot = makeFileRoot("out/testIterateGrid", outSuffix, D);
      openLogFile(outFileRoot + ".log");

      setupSystem(system, paramFileName);
      system.w().readBasis(wFileName);
      error = system.iterate();

      if (error) {
         Log::file() << "SCFT iterator failed to converge \n";
      } else {

         system.w().writeRGrid(outFileRoot + "_w.rf");
         system.c().writeRGrid(outFileRoot + "_c.rf");

         std::string refFileRoot;
         refFileRoot = makeFileRoot("ref/testIterate", outSuffix, D);
   
         FieldIo<D> const & fieldIo = system.domain().fieldIo();
         
         // Compare w fields
         std::string refFileBasis;
         std::string refFileRGrid = outFileRoot + "_ref_w.rf";
         if (convertRef) {
            refFileBasis = refFileRoot + "_w.bf";
            fieldIo.convertBasisToRGrid(refFileBasis, refFileRGrid);
         }
         wMaxDiff = readCompareWRGrid(system, refFileRGrid);
   
         // Optionally compare c fields
         if (compareCFields) {
            refFileBasis = refFileRoot + "_c.bf";
            refFileRGrid = outFileRoot + "_ref_c.rf";
            fieldIo.convertBasisToRGrid(refFileBasis, refFileRGrid);
            cMaxDiff = readCompareCRGrid(system, refFileRGrid);
         }

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
      UTIL_CHECK(system.scft().hasData());
      fHelmholtz = std::abs(fHelmholtz - system.scft().fHelmholtz());
      pressure   = std::abs(pressure - system.scft().pressure());
      if (verbose() > 0) {
         std::cout << "\nfHelmholtz diff = " << fHelmholtz;
         std::cout << "\npressure diff   = " << pressure;
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
      int error;
      testIterate(system,
                  "in/diblock/lam/param.rigid",
                   "in/diblock/lam/omega.ref",
                   "lam_rigid",
                   error,
                   wMaxDiff,
                   cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

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
      int error;
      testIterate(system,
                  "in/diblock/lam/param.flex",
                  "in/diblock/lam/omega.in",
                  "lam_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-7);

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

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref
   }

   void testIterate1D_lam_bead_flex()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/diblock/lam_bead/param.flex",
                  "in/diblock/lam_bead/w.bf",
                  "lam_bead_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 2.0E-8);
      TEST_ASSERT(cMaxDiff < 2.0E-8);

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

   void testIterate1D_lam_soln()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/solution/lam/param",
                  "in/solution/lam/w.bf",
                  "lam_soln",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
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

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to input
   }

   void testIterate1D_lam_open_soln()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      System<1> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/solution/lam_open/param",
                  "in/solution/lam_open/w.bf",
                  "lam_open_soln",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      //std::cout << "\n wMaxDiff = " << wMaxDiff;
      //std::cout << "\n cMaxDiff = " << cMaxDiff;
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-6);
      TEST_ASSERT(cMaxDiff < 1.0E-7);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = -2.75154924224e+01;
      double pressure   =  3.24415250699e+01;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

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
      int error;
      testIterate(system,
                  "in/blend/lam/param.open",
                  "in/blend/lam/w.bf",
                  "lam_open_blend",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = -3.30485085194e+00;
      double pressure = 5.10158598280e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to w.ref
   }

   #if 0
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
   #endif

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/diblock/hex/param.rigid",
                  "in/diblock/hex/omega.ref",
                  "hex_rigid",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz = 2.80222103795e+00;
      double pressure = 3.19716573940e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

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
      int error;
      testIterate(system,
                  "in/diblock/hex/param.flex",
                  "in/diblock/hex/omega.in",
                  "hex_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-7);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   2.80222103517e+00;
      double pressure =     3.19716584873e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.in as input, compared to omega.ref

      FieldIo<2> const & fieldIo = system.domain().fieldIo();
      fieldIo.scaleFieldsBasis("out/testIterateBasis2D_hex_flex_w.bf",
                               "out/testIterateBasis2D_hex_flex_w_scaled.bf",
                               0.01);
   }

   void testIterate2D_hex_bead_flex()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      double wMaxDiff, cMaxDiff;
      System<2> system;
      int error;
      testIterate(system,
                  "in/diblock/hex_bead/param.flex",
                  "in/diblock/hex_bead/w.bf",
                  "hex_bead_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 2.0E-7);
      TEST_ASSERT(cMaxDiff < 2.0E-7);

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

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/diblock/bcc/param.rigid",
                  "in/diblock/bcc/omega.ref",
                  "bcc_rigid",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(!error);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918380842e+00;
      double pressure =     4.03176984344e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used omega.ref as input, compared to input
   }

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      System<3> system;
      double wMaxDiff, cMaxDiff;
      int error;
      testIterate(system,
                  "in/diblock/bcc/param.flex",
                  "in/diblock/bcc/omega.in",
                  "bcc_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 2.0E-6);
      TEST_ASSERT(cMaxDiff < 3.0E-8);

      // Compare free energies to output of modified unit test from v1.1
      double fHelmholtz =   3.36918376624e+00;
      double pressure =     4.03176988267e+00;
      compareFreeEnergies(system, fHelmholtz, pressure);
      TEST_ASSERT(fHelmholtz < 1.0E-7);
      TEST_ASSERT(pressure < 1.0E-7);

      // Compare unit cell value
      double a = system.domain().unitCell().parameter(0);
      //std::cout << "a = " << Dbl(a,20,12) << std::endl;
      TEST_ASSERT(std::abs(a - 1.7593559883) < 1.0E-8);

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
      int error;
      testIterate(system,
                  "in/triblock/altGyr/param",
                  "in/triblock/altGyr/w.bf",
                  "altGyr_flex",
                  error,
                  wMaxDiff,
                  cMaxDiff);
      TEST_ASSERT(wMaxDiff < 1.0E-8);
      TEST_ASSERT(cMaxDiff < 1.0E-8);

      // Check stress value
      FSArray<double, 6> stress = computeStress(system);
      TEST_ASSERT(std::abs(stress[0]) < 1.0E-8);

      // v1.1 test used w.bf as input, compared to input

      // Output block concentrations - just tests that this doesn't crash
      std::string blockFile = "out/testIterate3D_altGyr_flex_block_c.rf";
      system.mixture().writeBlockCRGrid(blockFile);

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

};

TEST_BEGIN(AmIteratorGridTest)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_rigid)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_flex)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_bead_flex)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_soln)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_open_soln)
TEST_ADD(AmIteratorGridTest, testIterate1D_lam_open_blend)
//TEST_ADD(AmIteratorGridTest, testIterate1D_lam_open_shift)
TEST_ADD(AmIteratorGridTest, testIterate2D_hex_rigid)
TEST_ADD(AmIteratorGridTest, testIterate2D_hex_flex)
TEST_ADD(AmIteratorGridTest, testIterate2D_hex_bead_flex)
TEST_ADD(AmIteratorGridTest, testIterate3D_bcc_rigid)
TEST_ADD(AmIteratorGridTest, testIterate3D_bcc_flex)
TEST_ADD(AmIteratorGridTest, testIterate3D_altGyr_flex)
TEST_END(AmIteratorGridTest)

#endif
