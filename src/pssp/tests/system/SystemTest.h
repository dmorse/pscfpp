#ifndef PSSP_SYSTEM_TEST_H
#define PSSP_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/System.h>
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
      openInputFile("in/System1D", in);
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


    void testIterate3d()
   {   
     printMethod(TEST_FUNC);

     System<3> system;
     std::ifstream in; 
     openInputFile("in/System3D", in);
    
     //System<1> system;
     //std::ifstream in; 
     //openInputFile("in_check/System1D", in);
    
     system.readParam(in);
     in.close();
     std::ifstream command;
     openInputFile("in/Iterate3d", command);
     system.readCommands(command);
     command.close();

     int nMonomer = system.mixture().nMonomer();
     DArray<RField<3> > cFields_check;
     DArray<RField<3>  > cFields;
     cFields_check.allocate(nMonomer);
     cFields.allocate(nMonomer);
     double nx = system.mesh().size();
     for (int i = 0; i < nMonomer; ++i) {
         cFields_check[i].allocate(nx);
         cFields[i].allocate(nx);
     }   

     std::ifstream verify;
     openInputFile("contents/rho_rgrid_bcc", verify);

      std::string label;
      std::string uCell;
      std::string groupName;
      IntVec<3> nGrid;
      int nM; 
      int ver1, ver2;
      int dim;    
      int nCellParams;
      FArray<double,6> params;

      verify >> label;
      UTIL_CHECK(label == "format");
      verify >> ver1;
      verify >> ver2;

      verify >> label;
      UTIL_CHECK(label == "dim");
      verify >> dim;
      UTIL_CHECK(dim == 3);

      verify >> label;
      UTIL_CHECK(label == "crystal_system");
      verify >> uCell;

      verify >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify >> nCellParams;

      verify >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify >> params[i];
      }

      verify >> label;
      UTIL_CHECK(label == "group_name");
      verify >> groupName;

      verify >> label;
      UTIL_CHECK(label == "N_monomer");
      verify >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify >> label;
      UTIL_CHECK(label == "ngrid");
      verify >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());

      // Read Fields;
      MeshIterator<3> itr(system.mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify >>std::setprecision(15)>>cFields[i][itr.rank()];
         }
      }

     verify.close();

      std::ifstream verify_check;
      openInputFile("out/rho_check_bcc", verify_check);

      verify_check >> label;
      UTIL_CHECK(label == "format");
      verify_check >> ver1;
      verify_check >> ver2;

      verify_check >> label;
      UTIL_CHECK(label == "dim");
      verify_check >> dim;
      UTIL_CHECK(dim == 3);

      verify_check >> label;
      UTIL_CHECK(label == "crystal_system");
      verify_check >> uCell;

      verify_check >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify_check >> nCellParams;

      verify_check >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify_check >> params[i];
      }

      verify_check >> label;
      UTIL_CHECK(label == "group_name");
      verify_check >> groupName;

      verify_check >> label;
      UTIL_CHECK(label == "N_monomer");
      verify_check >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify_check >> label;
      UTIL_CHECK(label == "ngrid");
      verify_check >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());
 
      // Read Fields;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify_check >>std::setprecision(15)>>cFields_check[i][itr.rank()];
         }
      }

      verify_check.close();

      bool diff = true;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            if((std::abs(cFields_check[i][itr.rank()]-cFields[i][itr.rank()]))>1.0E-8){
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





    void testIterate1d()
   {
     printMethod(TEST_FUNC);

     System<1> system;
     std::ifstream in;
     openInputFile("in/System1D", in);
    
     system.readParam(in);
     in.close();
     std::ifstream command;
     openInputFile("in/Iterate1d", command);
     system.readCommands(command);
     command.close();

     int nMonomer = system.mixture().nMonomer();
     DArray<RField<1> > cFields_check;
     DArray<RField<1>  > cFields;
     cFields_check.allocate(nMonomer);
     cFields.allocate(nMonomer);
     double nx = system.mesh().size();
     for (int i = 0; i < nMonomer; ++i) {
         cFields_check[i].allocate(nx);
         cFields[i].allocate(nx);
     }

     std::ifstream verify;
     openInputFile("contents/rho_rgrid_lam", verify);

      std::string label;
      std::string uCell;
      std::string groupName;
      IntVec<1> nGrid;
      int nM;
      int ver1, ver2;
      int dim;      
      int nCellParams;
      FArray<double,6> params;

      verify >> label;
      UTIL_CHECK(label == "format");
      verify >> ver1;
      verify >> ver2;
 
      verify >> label;
      UTIL_CHECK(label == "dim");
      verify >> dim;
      UTIL_CHECK(dim == 1);

      verify >> label;
      UTIL_CHECK(label == "crystal_system");
      verify >> uCell;
      
      verify >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify >> nCellParams;

      verify >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify >> params[i];
      }
 
      verify >> label;
      UTIL_CHECK(label == "group_name");
      verify >> groupName;

      verify >> label;
      UTIL_CHECK(label == "N_monomer");
      verify >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify >> label;
      UTIL_CHECK(label == "ngrid");
      verify >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());

      // Read Fields;
      MeshIterator<1> itr(system.mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify >>std::setprecision(15)>>cFields[i][itr.rank()];
         }
      } 

     verify.close(); 


      std::ifstream verify_check;
      openInputFile("out/rho_check_lam", verify_check);

      verify_check >> label;
      UTIL_CHECK(label == "format");
      verify_check >> ver1;
      verify_check >> ver2;

      verify_check >> label;
      UTIL_CHECK(label == "dim");
      verify_check >> dim;
      UTIL_CHECK(dim == 1);

      verify_check >> label;
      UTIL_CHECK(label == "crystal_system");
      verify_check >> uCell;

      verify_check >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify_check >> nCellParams;

      verify_check >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify_check >> params[i];
      }   
 
      verify_check >> label;
      UTIL_CHECK(label == "group_name");
      verify_check >> groupName;

      verify_check >> label;
      UTIL_CHECK(label == "N_monomer");
      verify_check >> nM; 
      UTIL_CHECK(nM > 0); 
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify_check >> label;
      UTIL_CHECK(label == "ngrid");
      verify_check >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());
      
      // Read Fields;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify_check >>std::setprecision(15)>>cFields_check[i][itr.rank()];
         }
      }

      verify_check.close();

      bool diff = true;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            if((std::abs(cFields_check[i][itr.rank()]-cFields[i][itr.rank()]))>1.0E-8){
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




   void testIterate2d()
   {
     printMethod(TEST_FUNC);

     System<2> system;
     std::ifstream in;
     openInputFile("in/System2D", in);

     system.readParam(in);
     in.close();
     std::ifstream command;
     openInputFile("in/Iterate2d", command);
     system.readCommands(command);
     command.close();

     int nMonomer = system.mixture().nMonomer();
     DArray<RField<2> > cFields_check;
     DArray<RField<2>  > cFields;
     cFields_check.allocate(nMonomer);
     cFields.allocate(nMonomer);
     double nx = system.mesh().size();
     for (int i = 0; i < nMonomer; ++i) {
         cFields_check[i].allocate(nx);
         cFields[i].allocate(nx);
     }

     std::ifstream verify;
     openInputFile("contents/rho_rgrid_hex", verify);

      std::string label;
      std::string uCell;
      std::string groupName;
      IntVec<2> nGrid;
      int nM;
      int ver1, ver2;
      int dim;
      int nCellParams;
      FArray<double,6> params;

      verify >> label;
      UTIL_CHECK(label == "format");
      verify >> ver1;
      verify >> ver2;

      verify >> label;
      UTIL_CHECK(label == "dim");
      verify >> dim;
      UTIL_CHECK(dim == 2);

      verify >> label;
      UTIL_CHECK(label == "crystal_system");
      verify >> uCell;

      verify >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify >> nCellParams;

      verify >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify >> params[i];
      }
 
      verify >> label;
      UTIL_CHECK(label == "group_name");
      verify >> groupName;

      verify >> label;
      UTIL_CHECK(label == "N_monomer");
      verify >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify >> label;
      UTIL_CHECK(label == "ngrid");
      verify >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());

      // Read Fields;
      MeshIterator<2> itr(system.mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify >>std::setprecision(15)>>cFields[i][itr.rank()];
         }
      } 

     verify.close();


      std::ifstream verify_check;
      openInputFile("out/rho_check_hex", verify_check);

      verify_check >> label;
      UTIL_CHECK(label == "format");
      verify_check >> ver1;
      verify_check >> ver2;

      verify_check >> label;
      UTIL_CHECK(label == "dim");
      verify_check >> dim;
      UTIL_CHECK(dim == 2);

      verify_check >> label;
      UTIL_CHECK(label == "crystal_system");
      verify_check >> uCell;

      verify_check >> label;
      UTIL_CHECK(label == "N_cell_param");
      verify_check >> nCellParams;

      verify_check >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < nCellParams; ++i) {
         verify_check >> params[i];
      }
 
      verify_check >> label;
      UTIL_CHECK(label == "group_name");
      verify_check >> groupName;

      verify_check >> label;
      UTIL_CHECK(label == "N_monomer");
      verify_check >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == system.mixture().nMonomer());

      verify_check >> label;
      UTIL_CHECK(label == "ngrid");
      verify_check >> nGrid;
      UTIL_CHECK(nGrid == system.mesh().dimensions());

      // Read Fields;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
            verify_check >>std::setprecision(15)>>cFields_check[i][itr.rank()];
         }
      }

      verify_check.close();

      bool diff = true;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nM; ++i) {
           if((std::abs(cFields_check[i][itr.rank()]-cFields[i][itr.rank()]))>1.0E-8){
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




};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion3d)
TEST_ADD(SystemTest, testIterate3d)
TEST_ADD(SystemTest, testIterate1d)
TEST_ADD(SystemTest, testIterate2d)
TEST_END(SystemTest)

#endif
