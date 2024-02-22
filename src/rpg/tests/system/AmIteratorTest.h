#ifndef RPG_AM_ITERATOR_TEST_H
#define RPG_AM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/GpuResources.h>
#include <prdc/crystal/BFieldComparison.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class AmIteratorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

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
      copyFieldsBasis(b_wFields_check, system.w().basis());
    
     
      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate1D_lam_rigid_w.bf");
      system.writeCBasis("out/testIterate1D_lam_rigid_c.bf");

      // Compare result to original in basis format
      // DArray< DArray<double> > b_wFields;
      // copyFieldsBasis(b_wFields, system.w().basis());
      BFieldComparison comparison(1);
      // comparison.compare(b_wFields_check, b_wFields);
      comparison.compare(b_wFields_check, system.w().basis());
      if (verbose() > 0) {
         Log::file() << "\n";
         Log::file() << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);

   }

   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_flex.log");

      System<1> system;
      setupSystem<1>(system,"in/diblock/lam/param.flex"); 

      system.readWBasis("in/diblock/lam/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

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
      copyFieldsBasis(b_wFields, system.w().basis());

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
      openLogFile("out/testIterate1D_lam_soln.log");

      System<1> system;
      setupSystem<1>(system,"in/solution/lam/param"); 

      // Make reference copy of w fields
      system.readWBasis("in/solution/lam/w.bf");
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

      // iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_soln_w.bf");
      system.writeCBasis("out/testIterate1D_lam_soln_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

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

   void testIterate1D_lam_blend()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_blend.log");

      System<1> system;
      setupSystem<1>(system,"in/blend/lam/param.closed"); 

      // Make reference copy of w fields
      system.readWBasis("in/blend/lam/w.ref");
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/blend/lam/w.bf");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_blend_w.bf");
      system.writeCBasis("out/testIterate1D_lam_blend_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

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
      openLogFile("out/testIterate1D_lam_open_soln.log");

      System<1> system;
      setupSystem<1>(system,"in/solution/lam_open/param"); 

      // Make reference copy of w fields
      system.readWBasis("in/solution/lam_open/w.ref");
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

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
      copyFieldsBasis(b_wFields, system.w().basis());

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

   void testIterate1D_lam_open_blend()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_open_blend.log");

      System<1> system;
      setupSystem<1>(system,"in/blend/lam/param.open"); 

      // Make reference copy of w fields
      system.readWBasis("in/blend/lam/w.ref");
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

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
      copyFieldsBasis(b_wFields, system.w().basis());

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

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_rigid.log");

      System<2> system;
      setupSystem<2>(system,"in/diblock/hex/param.rigid"); 

      // Read reference solution
      system.readWBasis("in/diblock/hex/omega.ref");
      TEST_ASSERT(system.basis().isInitialized());

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

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
      copyFieldsBasis(b_wFields, system.w().basis());

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
      copyFieldsBasis(b_wFields_check, system.w().basis());

      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate2D_hex_flex_w.bf");
      system.writeCBasis("out/testIterate2D_hex_flex_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

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
      openLogFile("out/testIterate3D_bcc_rigid.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param.rigid"); 

      system.readWBasis("in/diblock/bcc/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_rigid_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_rigid_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

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

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_flex.log");

      System<3> system;
      setupSystem<3>(system,"in/diblock/bcc/param.flex"); 

      system.readWBasis("in/diblock/bcc/omega.ref");

      // Make reference copy of w fields
      DArray< DArray<double> > b_wFields_check;
      copyFieldsBasis(b_wFields_check, system.w().basis());

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_flex_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_flex_c.bf");

      // Get test result
      DArray< DArray<double> > b_wFields;
      copyFieldsBasis(b_wFields, system.w().basis());

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

TEST_BEGIN(AmIteratorTest)
TEST_ADD(AmIteratorTest, testIterate1D_lam_rigid)
TEST_ADD(AmIteratorTest, testIterate1D_lam_flex)
TEST_ADD(AmIteratorTest, testIterate1D_lam_soln)
TEST_ADD(AmIteratorTest, testIterate1D_lam_blend)
TEST_ADD(AmIteratorTest, testIterate1D_lam_open_blend)
TEST_ADD(AmIteratorTest, testIterate1D_lam_open_soln)
TEST_ADD(AmIteratorTest, testIterate2D_hex_rigid)
TEST_ADD(AmIteratorTest, testIterate2D_hex_flex)
TEST_ADD(AmIteratorTest, testIterate3D_bcc_rigid)
TEST_ADD(AmIteratorTest, testIterate3D_bcc_flex)

TEST_END(AmIteratorTest)

#endif
