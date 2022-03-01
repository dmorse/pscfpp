#ifndef PSPG_MIXTURE_TEST_H
#define PSPG_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/solvers/Mixture.h>
#include <pspg/solvers/Polymer.h>
#include <pspg/solvers/Block.h>
#include <pspg/solvers/Propagator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <util/math/Constants.h>

#include <pspg/math/GpuResources.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class MixtureTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Mixture<1> mixture;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      Mixture<1> mixture;

      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);
      in.close();
   }

   void testSolver1D()
   {

      printMethod(TEST_FUNC);
      Mixture<1> mixture;

      std::ifstream in;
      openInputFile("in/Mixture", in);
      mixture.readParam(in);
      
      UnitCell<1> unitCell;
      in >> unitCell;
      
      IntVec<1> d;
      in >> d;
      in.close();

      Mesh<1> mesh;
      mesh.setDimensions(d);
      mixture.setMesh(mesh);

      // Construct wavelist 
      WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeMinimumImages(mesh, unitCell);

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray<Mixture<1>::WField> d_wFields;
      DArray<Mixture<1>::CField> d_cFields;
      d_wFields.allocate(nMonomer);
      d_cFields.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         d_wFields[i].allocate(nx);
         d_cFields[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      cudaReal* wFields0 = new cudaReal[nx];
      cudaReal* wFields1 = new cudaReal[nx];
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields0[i] = 0.5 + cs;
         wFields1[i] = 0.5 - cs;
      }
      cudaMemcpy(d_wFields[0].cDField(), wFields0, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_wFields[1].cDField(), wFields1, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);

      mixture.compute(d_wFields, d_cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      #if 0
      std::cout << "Propagator(0,0), Q = " 
                << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      #endif 

   }

   void testSolver2D()
   {
      printMethod(TEST_FUNC);
      Mixture<2> mixture;

      std::ifstream in;
      openInputFile("in/Mixture2d", in);
      mixture.readParam(in);
      
      UnitCell<2> unitCell;
      in >> unitCell;
      
      IntVec<2> d;
      in >> d;
      in.close();

      Mesh<2> mesh;
      mesh.setDimensions(d);
      mixture.setMesh(mesh);

      // Construct wavelist 
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeMinimumImages(mesh, unitCell);

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray<Mixture<2>::WField> d_wFields;
      DArray<Mixture<2>::CField> d_cFields;
      d_wFields.allocate(nMonomer);
      d_cFields.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         d_wFields[i].allocate(nx);
         d_cFields[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      cudaReal* wFields0 = new cudaReal[nx];
      cudaReal* wFields1 = new cudaReal[nx];
      // Generate oscillatory wField
      int dx = mesh.dimension(0);
      int dy = mesh.dimension(1);
      double fx = 2.0*Constants::Pi/double(dx);
      double fy = 2.0*Constants::Pi/double(dy);
      double cx, cy;
      int k = 0;
      for (int i = 0; i < dx; ++i) {
         cx = cos(fx*double(i));
         for (int j = 0; j < dy; ++j) {
            cy = cos(fy*double(j));
            wFields0[k] = 0.5 + cx + cy;
            wFields1[k] = 0.5 - cx - cy;
            ++k;
         }
      }

      cudaMemcpy(d_wFields[0].cDField(), wFields0, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_wFields[1].cDField(), wFields1, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);

      mixture.compute(d_wFields, d_cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      #if 0
      std::cout << "Propagator(0,0), Q = " 
                << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      #endif 

   }

   void testSolver2D_hex()
   {
      printMethod(TEST_FUNC);
      Mixture<2> mixture;

      std::ifstream in;
      openInputFile("in/Mixture2d_hex", in);
      mixture.readParam(in);
      
      UnitCell<2> unitCell;
      in >> unitCell;
      
      IntVec<2> d;
      in >> d;
      in.close();

      Mesh<2> mesh;
      mesh.setDimensions(d);
      mixture.setMesh(mesh);

      // Construct wavelist 
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeMinimumImages(mesh, unitCell);

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray<Mixture<2>::WField> d_wFields;
      DArray<Mixture<2>::CField> d_cFields;
      d_wFields.allocate(nMonomer);
      d_cFields.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         d_wFields[i].allocate(nx);
         d_cFields[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      cudaReal* wFields0 = new cudaReal[nx];
      cudaReal* wFields1 = new cudaReal[nx];
      // Generate oscillatory wField
      int dx = mesh.dimension(0);
      int dy = mesh.dimension(1);
      double fx = 2.0*Constants::Pi/double(dx);
      double fy = 2.0*Constants::Pi/double(dy);
      double cx, cy;
      int k = 0;
      for (int i = 0; i < dx; ++i) {
         cx = cos(fx*double(i));
         for (int j = 0; j < dy; ++j) {
            cy = cos(fy*double(j));
            wFields0[k] = 0.5 + cx + cy;
            wFields1[k] = 0.5 - cx - cy;
            ++k;
         }
      }

      cudaMemcpy(d_wFields[0].cDField(), wFields0, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_wFields[1].cDField(), wFields1, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);

      mixture.compute(d_wFields, d_cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      #if 0
      std::cout << "Propagator(0,0), Q = " 
                << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      #endif 

   }

   void testSolver3D()
   {
      printMethod(TEST_FUNC);
      Mixture<3> mixture;

      std::ifstream in;
      openInputFile("in/Mixture3d", in);
      mixture.readParam(in);
      
      UnitCell<3> unitCell;
      in >> unitCell;
      
      IntVec<3> d;
      in >> d;
      in.close();

      Mesh<3> mesh;
      mesh.setDimensions(d);
      mixture.setMesh(mesh);

      // Construct wavelist 
      WaveList<3> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeMinimumImages(mesh, unitCell);

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray<Mixture<3>::WField> d_wFields;
      DArray<Mixture<3>::CField> d_cFields;
      d_wFields.allocate(nMonomer);
      d_cFields.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         d_wFields[i].allocate(nx);
         d_cFields[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      cudaReal* wFields0 = new cudaReal[nx];
      cudaReal* wFields1 = new cudaReal[nx];
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields0[i] = 0.5 + cs;
         wFields1[i] = 0.5 - cs;
      }

      cudaMemcpy(d_wFields[0].cDField(), wFields0, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_wFields[1].cDField(), wFields1, 
                  nx*sizeof(cudaReal), cudaMemcpyHostToDevice);

      mixture.compute(d_wFields, d_cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      #if 0
      std::cout << "Propagator(0,0), Q = " 
                << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << "Propagator(1,0), Q = " 
                << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << "Propagator(1,1), Q = " 
                << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << "Propagator(0,1), Q = " 
                << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      #endif 

   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor1D)
TEST_ADD(MixtureTest, testReadParameters1D)
TEST_ADD(MixtureTest, testSolver1D)
TEST_ADD(MixtureTest, testSolver2D)
TEST_ADD(MixtureTest, testSolver2D_hex)
TEST_ADD(MixtureTest, testSolver3D)
TEST_END(MixtureTest)

#endif
