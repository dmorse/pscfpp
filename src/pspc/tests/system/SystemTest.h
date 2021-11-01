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
      openInputFile("in_new/diblock/lam/param.flex", in);
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
      openInputFile("in_new/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w-fields (reference solution, solved by Fortran PSCF)
      system.readWBasis("in_new/diblock/lam/omega.in");

      // Copy w field components to wFields_check after reading
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();   

      // Round trip conversion basis -> rgrid -> basis, read result
      system.basisToRGrid("in_new/diblock/lam/omega.in",
                          "out/omega/conv/omega_rgrid_lam");
      system.rGridToBasis("out/omega/conv/omega_rgrid_lam",
                          "out/omega/conv/omega_conv_lam");
      system.readWBasis("out/omega/conv/omega_conv_lam");

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
      openInputFile("in_new/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in_new/diblock/hex/omega.in");

      // Store components in wFields_check for later comparison 
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Round trip basis -> rgrid -> basis, read resulting wField
      system.basisToRGrid("in_new/diblock/hex/omega.in",
                          "out/omega/conv/omega_rgrid_hex");
      system.rGridToBasis("out/omega/conv/omega_rgrid_hex",
                          "out/omega/conv/omega_conv_hex");
      system.readWBasis("out/omega/conv/omega_conv_hex");

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
      TEST_ASSERT(max < 1.0E-8);

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
      openInputFile("in_new/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      // Read w fields in system.wFields
      system.readWBasis("in_new/diblock/bcc/omega.in");
  
      // Store components of field as input 
      int nMonomer = system.mixture().nMonomer();
      int nStar = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Complete round trip basis -> rgrid -> basis
      system.basisToRGrid("in_new/diblock/bcc/omega.in",
                          "out/omega/conv/omega_rgrid_bcc");
      system.rGridToBasis("out/omega/conv/omega_rgrid_bcc",
                          "out/omega/conv/omega_conv_bcc");
      system.readWBasis("out/omega/conv/omega_conv_bcc");

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

   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_rigid.log"); 

      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in_new/diblock/lam/param.rigid", in); 
      system.readParam(in);
      in.close();

      // Read w fields
      system.readWBasis("in_new/diblock/lam/omega_lam");

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read w-fields, iterate and output solution
      system.iterate();
      system.writeWBasis("out/omega/domainOff/omega_lam");

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.07058e-08)){
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
      if (std::abs(system.mixture().stress(0) - 0.006583929) < 1.0E-8) {
         //0.006583929 is the stress calculated 
         //for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);

   }

   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_.log"); 
 
      System<1> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in; 
      openInputFile("in_new/diblock/lam/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in_new/diblock/lam/omega.in");
 
      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read input w-fields, iterate and output solution
      system.iterate();
      system.writeWBasis("out/omega/domainOn/omega_lam");

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)) {
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
      openInputFile("in_new/diblock/hex/param.rigid", in); 
      system.readParam(in);
      in.close();

      // Read reference solution
      system.readWBasis("in_new/diblock/hex/omega_hex");

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Read initial guess, iterate, output solution
      system.iterate();
      system.writeWBasis("out/omega/domainOff/omega_hex");

      // Compare current solution to reference solution
      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.60828e-07)) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.0e-07)) {
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
      if (std::abs(system.mixture().stress(0) - 0.010633960) < 1.0E-8) {
         // 0.010633960 is the stress calculated 
         // for this omega field for no stress relaxation using Fortran
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
      openInputFile("in_new/diblock/hex/param.flex", in);
      system.readParam(in);
      in.close();

      // Read reference solution (produced by Fortran code)
      system.readWBasis("in_new/diblock/hex/omega.in");

      // Save reference solution to wFields_check array
      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.iterate();
      system.writeWBasis("out/omega/domainOn/omega_hex");

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            // if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.58007e-07)) {
            if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.0e-07)) {
               // The above is the maximum error in the omega field.
               // Occurs for the first star
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
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
      openInputFile("in_new/diblock/bcc/param.rigid", in); 
      system.readParam(in);
      in.close();

      system.readWBasis("in_new/diblock/bcc/omega_bcc");

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.iterate();
      system.writeWBasis("out/omega/domainOff/omega_bcc");

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 1.02291e-07)) {
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.0e-07)) {
               // The above is the maximum error in the omega field.
               // Occurs for the second star.               
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
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
      if (std::abs(system.mixture().stress(0) - 0.005242863) < 1.0E-8) {
         //0.005242863 is the stress calculated for this omega field 
         //for no stress relaxation using Fortran
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
      openInputFile("in_new/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in_new/diblock/bcc/omega.in");

      int nMonomer = system.mixture().nMonomer();
      int ns = system.basis().nStar();
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      system.iterate();
      system.writeWBasis("out/omega/domainOn/omega_bcc");

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           //if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >=  1.09288e-07)) { 
           if ((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >=  5.0e-07)) { 
               // The above is the maximum error in the omega field.
               // Occurs for the second star.
               diff = false;
               std::cout <<"\n This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
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
TEST_ADD(SystemTest, testIterate1D_lam_rigid)
TEST_ADD(SystemTest, testIterate1D_lam_flex)
TEST_ADD(SystemTest, testIterate2D_hex_rigid)
TEST_ADD(SystemTest, testIterate2D_hex_flex)
TEST_ADD(SystemTest, testIterate3D_bcc_rigid)
TEST_ADD(SystemTest, testIterate3D_bcc_flex)

TEST_END(SystemTest)

#endif
