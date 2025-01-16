#ifndef RPG_MIXTURE_TEST_H
#define RPG_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Block.h>
#include <rpg/solvers/Propagator.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <util/math/Constants.h>

#include <pscf/cuda/GpuResources.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

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
      FFT<1> fft;
      fft.setup(d);
      mixture.setDiscretization(mesh, fft);

      // Construct wavelist 
      WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeAll();

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray< RField<1> > wFields;
      DArray< RField<1> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields_h[0][i] = 0.5 + cs;
         wFields_h[1][i] = 0.5 - cs;
      }
      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      if (verbose()) {
         std::cout << "Propagator(0,0), Q = " 
                   << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
         std::cout << "Propagator(1,0), Q = " 
                   << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
         std::cout << "Propagator(1,1), Q = " 
                   << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
         std::cout << "Propagator(0,1), Q = " 
                   << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      }

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
      FFT<2> fft;
      fft.setup(d);
      mixture.setDiscretization(mesh, fft);

      // Construct wavelist 
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeAll();

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray< RField<2> > wFields;
      DArray< RField<2> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
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
            wFields_h[0][k] = 0.5 + cx + cy;
            wFields_h[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }

      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      if (verbose()) {
         std::cout << "Propagator(0,0), Q = " 
                   << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
         std::cout << "Propagator(1,0), Q = " 
                   << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
         std::cout << "Propagator(1,1), Q = " 
                   << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
         std::cout << "Propagator(0,1), Q = " 
                   << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      }

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
      FFT<2> fft;
      fft.setup(d);
      mixture.setDiscretization(mesh, fft);

      // Construct wavelist 
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeAll();

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray< RField<2> > wFields;
      DArray< RField<2> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
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
            wFields_h[0][k] = 0.5 + cx + cy;
            wFields_h[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }

      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      if (verbose()) {
         std::cout << "Propagator(0,0), Q = " 
                   << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
         std::cout << "Propagator(1,0), Q = " 
                   << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
         std::cout << "Propagator(1,1), Q = " 
                   << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
         std::cout << "Propagator(0,1), Q = " 
                   << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      }

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
      FFT<3> fft;
      fft.setup(d);
      mixture.setDiscretization(mesh, fft);

      // Construct wavelist 
      WaveList<3> wavelist;
      wavelist.allocate(mesh, unitCell);
      wavelist.computeAll();

      // Setup unit cell
      mixture.setupUnitCell(unitCell, wavelist);

      int nMonomer = mixture.nMonomer();
      DArray< RField<3> > wFields;
      DArray< RField<3> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields_h[0][i] = 0.5 + cs;
         wFields_h[1][i] = 0.5 - cs;
      }

      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double Q = mixture.polymer(0).propagator(1, 0).computeQ();
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 1).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(0, 0).computeQ()));
      TEST_ASSERT(eq(Q, mixture.polymer(0).propagator(1, 1).computeQ()));

      if (verbose()) {
         std::cout << "Propagator(0,0), Q = " 
                   << mixture.polymer(0).propagator(0, 0).computeQ() << "\n";
         std::cout << "Propagator(1,0), Q = " 
                   << mixture.polymer(0).propagator(1, 0).computeQ() << "\n";
         std::cout << "Propagator(1,1), Q = " 
                   << mixture.polymer(0).propagator(1, 1).computeQ() << "\n";
         std::cout << "Propagator(0,1), Q = " 
                   << mixture.polymer(0).propagator(0, 1).computeQ() << "\n";
      }

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
