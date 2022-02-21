#ifndef PSPG_SYSTEM_TEST_H
#define PSPG_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/field/BFieldComparison.h>
#include <pspg/field/RDField.h>
#include <pspg/math/GpuResources.h>

//#include <pscf/mesh/MeshIterator.h>
//#include <util/format/Dbl.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class SystemTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  setVerbose(1); }

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
      system.setGpuResources(2, 16);

      std::ifstream in;
      openInputFile("in/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w-fields (reference solution, solved by Fortran PSCF)
      system.readWBasis("in/diblock/lam/omega.in");

      DArray< RDField<1> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      // Copy w field components to wFields_check after reading
      //DArray< RDField<1> > wFields_check;
      //wFields_check = system.wFields();

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in/diblock/lam/omega.in",
                          "out/testConversion1D_lam_w.rf");
      system.rGridToBasis("out/testConversion1D_lam_w.rf",
                          "out/testConversion1D_lam_w.bf");
      system.readWBasis("out/testConversion1D_lam_w.bf");


      DArray< RDField<1> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }   

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     } 

      // Compare result to original
      BFieldComparison comparison;
      comparison.compare(wFields_check, wFields_test, nStar);
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
      system.setGpuResources(32, 32);

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in/diblock/hex/omega.in");

      DArray< RDField<2> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }   

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }   

      // Store components in wFields_check for later comparison
      //DArray< RDField<2>  > wFields_check;
      //wFields_check = system.wFields();

      // Round trip basis -> rgrid -> basis, read resulting wField
      system.basisToRGrid("in/diblock/hex/omega.in",
                          "out/testConversion2D_hex_w.rf");

      system.rGridToBasis("out/testConversion2D_hex_w.rf",
                          "out/testConversion2D_hex_w.bf");
      system.readWBasis("out/testConversion2D_hex_w.bf");

      // Check symmetry of rgrid representation
      //bool hasSymmetry
      // = system.checkRGridFieldSymmetry("out/testConversion2D_hex_w.rf");
      //TEST_ASSERT(hasSymmetry);

      DArray< RDField<2> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }   

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison;
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
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
      system.setGpuResources(32, 1024);

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();
      // Read w fields in system.wFields
      system.readWBasis("in/diblock/bcc/omega.in");

      DArray< RDField<3> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      // Store components of field as input
      //DArray< RDField<3> > wFields_check;
      //wFields_check = system.wFields();

      // Complete round trip basis -> rgrid -> basis
      system.basisToRGrid("in/diblock/bcc/omega.in",
                          "out/testConversion3D_bcc_w.rf");
      system.rGridToBasis("out/testConversion3D_bcc_w.rf",
                          "out/testConversion3D_bcc_w.bf");
      system.readWBasis("out/testConversion3D_bcc_w.bf");

      // Check symmetry of rgrid representation
      // bool hasSymmetry
      // = system.checkRGridFieldSymmetry("out/testConversion3D_bcc_w.rf");
      //TEST_ASSERT(hasSymmetry);

      // Compare result to original
      DArray< RDField<3> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
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
*      system.setGpuResources(32, 1024);
*
*      // Read system parameter file
*      std::ifstream in;
*      openInputFile("in/diblock/bcc/param.flex", in);
*      system.readParam(in);
*      in.close();
*
*      system.readWBasis("in/diblock/bcc/omega.in");
*      bool hasSymmetry = system.fieldIo().hasSymmetry(system.wFieldRGrid(0));
*      TEST_ASSERT(hasSymmetry);
*
*      // Intentionally mess up the field, check that symmetry is destroyed
*      system.wFieldRGrid(0)[23] += 0.1;
*      hasSymmetry = system.fieldIo().hasSymmetry(system.wFieldRGrid(0));
*      TEST_ASSERT(!hasSymmetry);
*
*   }
*/

   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate1D_lam_rigid.log");
      system.setGpuResources(1, 32);

      std::ifstream in;
      openInputFile("in/diblock/lam/param.rigid", in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in/diblock/lam/omega.ref");

      DArray< RDField<1> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }
     
      system.readWBasis("in/diblock/lam/omega.in");
      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate1D_lam_rigid_w.bf");
      system.writeCBasis("out/testIterate1D_lam_rigid_c.bf");

      DArray< RDField<1> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
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

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate1D_lam_flex.log");
      system.setGpuResources(1, 32);

      std::ifstream in;
      openInputFile("in/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/lam/omega.ref");

      DArray< RDField<1> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      //DArray< RDField<1> > wFields_check;
      //wFields_check = system.wFields();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate1D_lam_flex_w.bf");
      system.writeCBasis("out/testIterate1D_lam_flex_c.bf");

      DArray< RDField<1> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

   }

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate2D_hex_rigid.log");
      system.setGpuResources(8, 128);

      std::ifstream in;
      openInputFile("in/diblock/hex/param.rigid", in);
      system.readParam(in);
      in.close();

      // Read reference solution
      system.readWBasis("in/diblock/hex/omega.ref");

      DArray< RDField<2> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      //DArray< RDField<2> > wFields_check;
      //wFields_check = system.wFields();

      // Read initial guess, iterate, output solution
      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate2D_hex_rigid_w.bf");
      system.writeCBasis("out/testIterate2D_hex_rigid_c.bf");

      DArray< RDField<2> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

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

      System<2> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate2D_hex_flex.log");
      system.setGpuResources(8, 128);

      // Read parameter file
      std::ifstream in;
      openInputFile("in/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read reference solution (produced by Fortran code)
      system.readWBasis("in/diblock/hex/omega.ref");

      DArray< RDField<2> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      // Save reference solution to wFields_check array
      //DArray< RDField<2> > wFields_check;
      //wFields_check = system.wFields();

      system.readWBasis("in/diblock/hex/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate2D_hex_flex_w.bf");
      system.writeCBasis("out/testIterate2D_hex_flex_c.bf");

      DArray< RDField<2> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate3D_bcc_rigid.log");
      system.setGpuResources(128, 256);

      std::ifstream in;
      openInputFile("in/diblock/bcc/param.rigid", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/bcc/omega.ref");

      DArray< RDField<3> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      //DArray< RDField<3> > wFields_check;
      //wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_rigid_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_rigid_c.bf");

      DArray< RDField<3> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 7E-7);

      

   }

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testIterate3D_bcc_flex.log");
      system.setGpuResources(128, 256);

      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/bcc/omega.ref");

      DArray< RDField<3> > d_wFields_check;
      d_wFields_check = system.wFields();

      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray<cudaReal*> wFields_check;
      wFields_check.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_check[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_check[i], d_wFields_check[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      //DArray< RDField<3> > wFields_check;
      //wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      };
      system.writeWBasis("out/testIterate3D_bcc_flex_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_flex_c.bf");

      DArray< RDField<3> > wFieldsGpu_test;
      wFieldsGpu_test = system.wFields();

      DArray<cudaReal*> wFields_test;
      wFields_test.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         wFields_test[i] = new cudaReal[nStar];
      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(wFields_test[i], wFieldsGpu_test[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }


      // Compare result to original
      BFieldComparison comparison (1);
      comparison.compare(wFields_check, wFields_test, nStar);
      if (verbose()>0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);


   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
// TEST_ADD(SystemTest, testCheckSymmetry3D_bcc)
TEST_ADD(SystemTest, testIterate1D_lam_rigid)
TEST_ADD(SystemTest, testIterate1D_lam_flex)
// TEST_ADD(SystemTest, testIterate2D_hex_rigid)
// TEST_ADD(SystemTest, testIterate2D_hex_flex)
TEST_ADD(SystemTest, testIterate3D_bcc_rigid)
TEST_ADD(SystemTest, testIterate3D_bcc_flex)


TEST_END(SystemTest)

#endif
