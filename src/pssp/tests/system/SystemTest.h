#ifndef PSSP_SYSTEM_TEST_H
#define PSSP_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/System.h>

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
      openInputFile("in/System", in);
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
      printMethod(TEST_FUNC);
      System<3> system;

      std::ifstream in;
      openInputFile("in/System3D", in);
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
    // System<3> system;
    // std::ifstream in;
    // openInputFile("in/System3D", in);
    
     System<1> system;
     std::ifstream in; 
     openInputFile("in/System1D", in);
     
     system.readParam(in);
     in.close();
     std::ifstream command;
     openInputFile("in/Iterate", command);
     system.readCommands(command);
     command.close();

   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion3d)
TEST_ADD(SystemTest, testIterate3d)
TEST_END(SystemTest)

#endif
