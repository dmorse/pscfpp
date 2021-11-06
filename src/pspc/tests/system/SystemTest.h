#ifndef PSPC_SYSTEM_TEST_H
#define PSPC_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pscf/mesh/MeshIterator.h>

//#include <pspc/iterator/AmIterator.h>
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
   {}

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
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();   

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in/diblock/lam/omega.in",
                          "out/diblock/lam/omega.rgrid");
      system.rGridToBasis("out/diblock/lam/omega.rgrid",
                          "out/diblock/lam/omega.conv");
      system.readWBasis("out/diblock/lam/omega.conv");

      // Compare result to original
      double err;
      double max = 0.0;
      std::cout << std::endl;   
      for (int j = 0; j < nStar; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            err = wFields_check[i][j] - system.wFields()[i][j];
            err = std::abs(err);
            //std::cout << Dbl(wFields_check[i][j],15,8)  << "  ";
            //std::cout << Dbl(system.wFields()[i][j],15,8) << "  ";
            if (err > max) {
               max = err;
            }
         }
         //std::cout << std::endl;   
      }
      std::cout << "Max error = " << max << std::endl;  
      TEST_ASSERT(max < 1.0E-8);
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
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Round trip basis -> rgrid -> basis, read resulting wField
      system.basisToRGrid("in/diblock/hex/omega.in",
                          "out/diblock/hex/omega.rgrid");
      system.rGridToBasis("out/diblock/hex/omega.rgrid",
                          "out/diblock/hex/omega.conv");
      system.readWBasis("out/diblock/hex/omega.conv");

      // Check symmetry of rgrid representation
      bool hasSymmetry
       = system.checkRGridFieldSymmetry("out/diblock/hex/omega.rgrid");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      double err;
      double max = 0.0;
      std::cout << std::endl;   
      for (int j = 0; j < nStar; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            err = wFields_check[i][j] - system.wFields()[i][j];
            err = std::abs(err);
            //std::cout << wFields_check[i][j] 
            //          << "  " << system.wFields()[i][j];
            if (err > max) {
               max = err;
            }
         }
         //std::cout << std::endl;   
      }   
      std::cout << "Max error = " << max << std::endl;  
      TEST_ASSERT(max < 1.0E-9);

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
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Complete round trip basis -> rgrid -> basis
      system.basisToRGrid("in/diblock/bcc/omega.in",
                          "out/diblock/bcc/omega.rgrid");
      system.rGridToBasis("out/diblock/bcc/omega.rgrid",
                          "out/diblock/bcc/omega.conv");
      system.readWBasis("out/diblock/bcc/omega.conv");

      // Check symmetry of rgrid representation
      bool hasSymmetry
       = system.checkRGridFieldSymmetry("out/diblock/bcc/omega.rgrid");
      TEST_ASSERT(hasSymmetry);

      // Compare result to original
      double err;
      double max = 0.0;
      std::cout << std::endl;   
      for (int j = 0; j < nStar; ++j) {
         //std::cout << j << "  ";
         for (int i = 0; i < nMonomer; ++i) {
            err = wFields_check[i][j] - system.wFields()[i][j];
            err = std::abs(err);
            //std::cout << Dbl(wFields_check[i][j],15,8)  << "  ";
            //std::cout << Dbl(system.wFields()[i][j],15,8) << "  ";
            if (err > max) {
               max = err;
            }
         }
         //std::cout << std::endl;   
      }
      std::cout << "Max error = " << max << std::endl;  
      TEST_ASSERT(max < 1.0E-8);

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
      std::cout<<"Read Parameters";
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in/diblock/lam/omega.ref");

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/lam/omega.rigid");
      system.writeCBasis("out/diblock/lam/rho.rigid");

      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 1e-10)){
                // The above is the minimum error in the omega field.
                // Occurs for the first star                 
                diff = false;
                std::cout <<"\n This is error for break:"<<
                   (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
                std::cout << "star index = " << j << std::endl;
                break;
            } else {
                diff = true;
            }
         }    
         if (diff == false) {
            break;
         }    
      }    
      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-8) {
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);

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
 
      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read input w-fields, iterate and output solution
      system.readWBasis("in/diblock/lam/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/lam/omega.flex");
      system.writeCBasis("out/diblock/lam/rho.flex");

      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-10)) {
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"star index = "<< j << std::endl;
               break;
            }   
            else
               diff = true;
         }   
         if (diff==false) {
            break;
         }
      }   
      TEST_ASSERT(diff);
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

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read initial guess, iterate, output solution
      system.readWBasis("in/diblock/hex/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/hex/omega.rigid");
      system.writeCBasis("out/diblock/hex/rho.rigid");

      // Compare current solution to reference solution
      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.60828e-07)) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 3.0E-7)) {
               // The above is the minimum error in the omega field.
               // Occurs for the first star            
               diff = false;
               std::cout << "\n This is error for break:" << 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"star index = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }    
         if (diff==false) {
            break;
         }    
      }    
      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-8) {
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
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
      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/hex/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/hex/omega.flex");
      system.writeCBasis("out/diblock/hex/rho.flex");

      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            // if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.58007e-07)) {
            if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 3.0e-7)) {
               // The above is the maximum error in the omega field.
               // Occurs for the first star
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"star index = "<< j << std::endl;
               break;
            } else {
               diff = true;
            }
         }
         if (diff==false) {
            break;
         }
      }
      TEST_ASSERT(diff);
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

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/bcc/omega.rigid");
      system.writeCBasis("out/diblock/bcc/rho.rigid");

      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 1.02291e-07)) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.0e-07)) {
               // The above is the maximum error in the omega field.
               // Occurs for the second star.               
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"star index = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }
         if (diff==false) {
            break;
         }
      }
      bool stress = false;
      if (std::abs(system.mixture().stress(0)) < 1.0E-7) {
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
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

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.readWBasis("in/diblock/bcc/omega.in");
      system.iterate();
      system.writeWBasis("out/diblock/bcc/omega.flex");
      system.writeCBasis("out/diblock/bcc/rho.flex");

      bool diff = true;
      for (int j = 1; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >=  1.09288e-07)) { 
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >=  5.0e-07)) { 
               // The above is the maximum error in the omega field.
               // Occurs for the second star.
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"star index = "<< j << std::endl;
               break;
            }
            else
               diff = true;
         }
         if (diff==false) {

            break;
         }
      }
      TEST_ASSERT(diff);
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
