#ifndef RPG_SYSTEM_TEST_H
#define RPG_SYSTEM_TEST_H

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
      setupSystem<1>(system,"in/diblock/lam/param.flex"); 
   }

   void testConversion1D_lam()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testConversion1D_lam.log");
      
      System<1> system;
      setupSystem<1>(system,"in/diblock/lam/param.flex"); 

      // Read w-fields (reference solution, solved by Fortran PSCF)
      system.readWBasis("in/diblock/lam/omega.in");

      // Get reference field
      DArray< DArray<double> > b_wFields_check;
      //copyFieldsBasis(b_wFields_check, system.w().basis());
      b_wFields_check = system.w().basis();

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in/diblock/lam/omega.in",
                          "out/testConversion1D_lam_w.rf");
      system.rGridToBasis("out/testConversion1D_lam_w.rf",
                          "out/testConversion1D_lam_w.bf");
      system.readWBasis("out/testConversion1D_lam_w.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      //copyFieldsBasis(b_wFields, system.w().basis());
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);
      if (verbose()>0) {
         Log::file() << "\n";
         Log::file() << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testConversion2D_hex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testConversion2D_hex.log");

      System<2> system;
      setupSystem<2>(system,"in/diblock/hex/param.flex"); 

      // Read w fields
      system.readWBasis("in/diblock/hex/omega.in");

      // Get reference field
      DArray< DArray<double> > b_wFields_check;
      //copyFieldsBasis(b_wFields_check, system.w().basis());
      b_wFields_check = system.w().basis();

      // Round trip basis -> rgrid -> basis, read resulting wField
      system.basisToRGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.rf");

      system.rGridToBasis("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w.bf");
      system.readWBasis("out/testConversion2D_hex_w.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      //copyFieldsBasis(b_wFields, system.w().basis());
      b_wFields = system.w().basis();

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);
      if (verbose()>0) {
         Log::file() << "\n";
         Log::file() << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testConversion3D_bcc()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testConversion3D_bcc.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param.flex"); 


      // Read w fields in system.wFields
      system.readWBasis("in/diblock/bcc/omega.in");

      // Get reference field
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

      // Complete round trip basis -> rgrid -> basis
      system.basisToRGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.rf");
      system.rGridToBasis("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w.bf");
      system.readWBasis("out/testConversion3D_bcc_w.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(b_wFields_check, b_wFields);
      if (verbose()>0) {
         Log::file() << "\n";
         Log::file() << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

/*   void testCheckSymmetry3D_bcc()
*   {
*      printMethod(TEST_FUNC);
*      System<3> system;
*      system.fileMaster().setInputPrefix(filePrefix());
*      system.fileMaster().setOutputPrefix(filePrefix());
*
*      openLogFile("out/testSymmetry3D_bcc.log");
*
*      // Read system parameter file
*      std::ifstream in;
*      openInputFile("in/diblock/bcc/param.flex", in);
*      system.readParam(in);
*      in.close();
*
*      system.readWBasis("in/diblock/bcc/omega.in");
*      bool hasSymmetry = system.fieldIo().hasSymmetry(system.w().rgrid(0));
*      TEST_ASSERT(hasSymmetry);
*
*      // Intentionally mess up the field, check that symmetry is destroyed
*      system.w().rgrid(0)[23] += 0.1;
*      hasSymmetry = system.fieldIo().hasSymmetry(system.w().rgrid(0));
*      TEST_ASSERT(!hasSymmetry);
*
*   }
*/

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

   void copyFieldsBasis(DArray< DArray<double> > & out, 
                        DArray< DArray<double> > const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int nField = in.capacity();
      UTIL_CHECK(nField > 0);

      // If array out is not allocated, allocate
      if (!out.isAllocated()) {
         out.allocate(nField);
      }
      int nPoint = in[0].capacity();
      for (int i = 0; i < nField; i++) {
         UTIL_CHECK(in[i].capacity() == nPoint);
         if (!out[i].isAllocated()) {
            out[i].allocate(nPoint);
         } else {
            UTIL_CHECK(out[i].capacity() == nPoint);
         }
      }

      // Copy arrays
      for (int i = 0; i < nField; i++) {
         UTIL_CHECK(out[i].capacity() == in[i].capacity());
         out[i] = in[i];
      }
   }

   #if 0
   template <int D>
   void copyFieldsRGrid(DArray< RField<D> > & out, 
                        DArray< RField<D> > const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int nField = in.capacity();
      UTIL_CHECK(nField > 0);

      // If array out is not allocated, allocate
      if (!out.isAllocated()) {
         out.allocate(nField);
      }
      IntVec<D> dim = in[0].meshDimensions();
      for (int i = 0; i < nField; i++) {
         UTIL_CHECK(in[i].meshDimensions() == dim);
         if (!out[i].isAllocated()) {
            out[i].allocate(dim);
         } else {
            UTIL_CHECK(out[i].meshDimensions() == dim);
         }
      }

      // Copy arrays
      for (int i = 0; i < nField; i++) {
         UTIL_CHECK(out[i].capacity() == in[i].capacity());
         UTIL_CHECK(out[i].meshDimensions() == in[i].meshDimensions());
         out[i] = in[i];
      }

   }
   #endif

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
//// TEST_ADD(SystemTest, testCheckSymmetry3D_bcc)

TEST_END(SystemTest)

#endif
