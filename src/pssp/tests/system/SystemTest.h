#ifndef PSSP_SYSTEM_TEST_H
#define PSSP_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/System.h>
#include <pssp/iterator/AmIterator.h>
#include <pscf/mesh/MeshIterator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class SystemTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;

      std::ifstream in;
      openInputFile("in/domainOn/System1D", in);
      system.readParam(in);
      in.close();
   }
 
   #if 0
   void testSolver1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;

      std::ifstream in;
      openInputFile("in_check/System", in);
      system.readParam(in);
      UnitCell<1> unitCell;
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      in.close();

      Mesh<1> mesh;
      mesh.setDimensions(d);
      system.setMesh(mesh);
      system.setupUnitCell(unitCell);

      std::cout << "\n";
      system.writeParam(std::cout);
      std::cout << "unitCell  " << unitCell << std::endl;
      std::cout << "mesh      " << mesh.dimensions() << std::endl;

      int nMonomer = system.nMonomer();
      DArray<System<1>::WField> wFields;
      DArray<System<1>::CField> cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(nx);
         cFields[i].allocate(nx);
      }

      double cs;
      for (int i = 0; i < nx; ++i) {
         //cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         //cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }

      system.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      std::cout << "Propagator(0,0), Q = " 
                << system.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << system.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << system.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << system.polymer(0).propagator(0, 1).computeQ() << "\n";

      #if 0
      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      #endif
      
   }
   #endif

   void testConversion3d()
   {
      /*printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in;
      openInputFile("in_check/System3D", in);

      System<1> system; 
      std::ifstream in;
      openInputFile("in_check/System1D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in_check/Conversion", command);
      system.readCommands(command);
      command.close();*/

      printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in; 
      openInputFile("in/System3D", in);

     // System<1> system; 
     // std::ifstream in; 
     // openInputFile("in/System1D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/Conversion", command);
      system.readCommands(command);
      command.close();

   }


   void testConversion_BCC() 
   {   
      printMethod(TEST_FUNC);
 
      System<3> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System3D", in);
          
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_3d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){         
            wFields_check[i][j] = system.wFields() [i] [j];
         }
      }

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_3d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
               break;
            }
            else
               diff = true;
         }
         if(diff==false)
            break;
      }
      TEST_ASSERT(diff);
   }   


   void testConversion_hex() 
   {   
      printMethod(TEST_FUNC);
 
      System<2> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System2D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_2d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_2d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
               break;
            }   
            else
               diff = true;
         }   
         if(diff==false)
            break;
      }   
      TEST_ASSERT(diff);
   }  


   void testConversion_lam() 
   {   
      printMethod(TEST_FUNC);
 
      System<1> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System1D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_1d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_1d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
               break;
            }   
            else
               diff = true;
         }   
         if(diff==false)
            break;
      }   
      TEST_ASSERT(diff);
   }  

   void testIterate_bcc_domainOn()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System3D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_bcc", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate3d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }
            else
               diff = true;
         }
         if(diff==false){

            break;
         }
      }
      TEST_ASSERT(diff);

   }


  void testIterate_bcc_domainOff()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in;
      openInputFile("in/domainOff/System3D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_bcc", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate3d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }
         if(diff==false){
            break;
         }
      }
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.005242863) < 1.0E-8){
         //0.005242863 is the stress calculated for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }

   void testIterate_hex_domainOn()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      std::ifstream in;
      openInputFile("in/domainOn/System2D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_hex", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j];
         }
      }

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate2d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }
            else
               diff = true;
         }
         if(diff==false){

            break;
         }
      }
      TEST_ASSERT(diff);
   }

   void testIterate_hex_domainOff()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      std::ifstream in;
      openInputFile("in/domainOff/System2D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_hex", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate2d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }    
         if(diff==false){
            break;
         }    
      }    
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.010633960) < 1.0E-8){
         //0.010633960 is the stress calculated for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }


   void testIterate_lam_domainOn()
   {
      printMethod(TEST_FUNC);
 
      System<1> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System1D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_lam", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate1d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }   
            else
               diff = true;
         }   
         if(diff==false){
            break;
         }
      }   
      TEST_ASSERT(diff);
   }

   void testIterate_lam_domainOff()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      std::ifstream in;
      openInputFile("in/domainOff/System1D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_lam", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate1d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
            std::cout <<"This is error for break:"<< (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
            std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }    
         if(diff==false){
            break;
         }    
      }    
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.006583929) < 1.0E-8){
         //0.006583929 is the stress calculated for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }

   void testIterate_BCC() 
   {   
      printMethod(TEST_FUNC);
 
      System<3> system;

      // Read parameter file
      std::ifstream in; 
      openInputFile("in/bcc/param", in);
      system.readParam(in);
      in.close();

      openInputFile("in/bcc/omega", in);
      system.readFields(in, system.wFields());
      in.close();

      int err = system.iterator().solve();
      TEST_ASSERT(!err);

   }


};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)

TEST_ADD(SystemTest, testConversion_BCC)
TEST_ADD(SystemTest, testConversion_hex)
TEST_ADD(SystemTest, testConversion_lam)
TEST_ADD(SystemTest, testIterate_bcc_domainOn)
TEST_ADD(SystemTest, testIterate_hex_domainOn)
TEST_ADD(SystemTest, testIterate_lam_domainOn)
TEST_ADD(SystemTest, testIterate_bcc_domainOff)
TEST_ADD(SystemTest, testIterate_hex_domainOff)
TEST_ADD(SystemTest, testIterate_lam_domainOff)

//TEST_ADD(SystemTest, testIterate_BCC)
TEST_END(SystemTest)

#endif
