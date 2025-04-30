#ifndef PRDC_CPU_WAVE_LIST_TEST_H
#define PRDC_CPU_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/WaveList.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;

class CpuWaveListTest : public UnitTest
{

private:

   // Error tolerance for array equality
   constexpr static double tolerance_ = 1E-10;

   Mesh<1> mesh1;
   Mesh<2> mesh2;
   Mesh<3> mesh3;

   UnitCell<1> cell1;
   UnitCell<2> cell2;
   UnitCell<3> cell3;

   IntVec<1> kMeshDims1;
   IntVec<2> kMeshDims2;
   IntVec<3> kMeshDims3;

   int kSize1, kSize2, kSize3;

public:

   void setUp()
   {
      // set up 1D mesh
      IntVec<1> meshDims1; 
      meshDims1[0] = 32; 
      mesh1.setDimensions(meshDims1);
      IntVec<1> kMeshDims1;
      kMeshDims1[0] = (meshDims1[0] / 2) + 1;
      kSize1 = kMeshDims1[0];

      // Set up 2D meshes
      IntVec<2> meshDims2; 
      meshDims2[0] = 32; 
      meshDims2[1] = 48; 
      mesh2.setDimensions(meshDims2);
      kMeshDims2 = meshDims2;
      kMeshDims2[1] = meshDims2[1] / 2 + 1;
      kSize2 = kMeshDims2[0] * kMeshDims2[1];

      // Set up 3D meshes
      IntVec<3> meshDims3; 
      meshDims3[0] = 24; 
      meshDims3[1] = 64; 
      meshDims3[2] = 21;
      mesh3.setDimensions(meshDims3);
      kMeshDims3 = meshDims3;
      kMeshDims3[2] = meshDims3[2] / 2 + 1;
      kSize3 = kMeshDims3[0] * kMeshDims3[1] * kMeshDims3[2];

      // set up 1D, 2D, and 3D unit cell objects
      std::ifstream in1, in2, in3;
      openInputFile("in/Lamellar", in1);
      in1 >> cell1;

      openInputFile("in/Oblique", in2);
      in2 >> cell2;

      openInputFile("in/Triclinic", in3);
      in3 >> cell3;
   }

   void tearDown()
   {}

   void testAllocate()
   {
      printMethod(TEST_FUNC);

      WaveList<1> wavelist1;
      wavelist1.allocate(mesh1, cell1);
      TEST_ASSERT(wavelist1.isAllocated());
      TEST_ASSERT(!wavelist1.hasMinimumImages());

      WaveList<2> wavelist2;
      wavelist2.allocate(mesh2, cell2);
      TEST_ASSERT(wavelist2.isAllocated());
      TEST_ASSERT(!wavelist2.hasMinimumImages());

      WaveList<3> wavelist3;
      wavelist3.allocate(mesh3, cell3);
      TEST_ASSERT(wavelist3.isAllocated());
      TEST_ASSERT(!wavelist3.hasMinimumImages());
   }

   void testComputeMinimumImages1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] = 17);
      TEST_ASSERT(kSize1 = 17);

      // set up wavelist object
      WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute minimum images (and ksq) on device, transfer to host
      wavelist.computeMinimumImages(); 
      DArray< IntVec<1> > const & minImages_h = wavelist.minImages();
      RField<1> const & ksq_h = wavelist.kSq();

      // Compute minimum images (and ksq) on host and compare
      IntVec<1> temp, vec;
      MeshIterator<1> iter;
      double ksq;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         ksq = cell1.ksq(vec);
         
         TEST_ASSERT(vec == minImages_h[iter.rank()]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }

   }

   void testComputeMinimumImages2D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh2.dimension(0) == 32);
      TEST_ASSERT(mesh2.dimension(1) == 48);
      TEST_ASSERT(kMeshDims2[0] == 32);
      TEST_ASSERT(kMeshDims2[1] == 25);
      TEST_ASSERT(kSize2 == 32*25);

      // Set up wavelist object
      WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute minimum images (and ksq) on device, transfer to host
      wavelist.computeMinimumImages(); 
      DArray< IntVec<2> > const & minImages_h = wavelist.minImages();
      RField<2> ksq_h = wavelist.kSq();
      TEST_ASSERT(minImages_h.capacity() == kSize2);
      TEST_ASSERT(ksq_h.capacity() == kSize2);

      // compute minimum images (and ksq) on host and compare
      IntVec<2> temp, vec;
      MeshIterator<2> iter;
      double ksq;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
         ksq = cell2.ksq(vec);
         
         TEST_ASSERT(vec == minImages_h[iter.rank()]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputeMinimumImages3D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh3.dimension(0) == 24);
      TEST_ASSERT(mesh3.dimension(1) == 64);
      TEST_ASSERT(mesh3.dimension(2) == 21);
      TEST_ASSERT(kMeshDims3[0] == 24);
      TEST_ASSERT(kMeshDims3[1] == 64);
      TEST_ASSERT(kMeshDims3[2] == 11);
      TEST_ASSERT(kSize3 == 24*64*11);

      // Set up wavelist object
      WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);

      // Compute minimum images (and ksq) on device, transfer to host
      wavelist.computeMinimumImages(); 
      DArray< IntVec<3> > const & minImages_h = wavelist.minImages();
      RField<3> const & ksq_h = wavelist.kSq();

      // Compute minimum images (and ksq) on host and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double ksq;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         ksq = cell3.ksq(vec);

         TEST_ASSERT(vec == minImages_h[iter.rank()]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] = 17);
      TEST_ASSERT(kSize1 = 17);

      // set up wavelist object
      WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<1> const & ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq 
      RField<1> const & ksq_h2 = wavelist.kSq();

      // Compute kSq on host and compare
      IntVec<1> temp, vec;
      double ksq;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         ksq = cell1.ksq(vec);

         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(ksq - ksq_h2[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq2D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh2.dimension(0) == 32);
      TEST_ASSERT(mesh2.dimension(1) == 48);
      TEST_ASSERT(kMeshDims2[0] == 32);
      TEST_ASSERT(kMeshDims2[1] == 25);
      TEST_ASSERT(kSize2 == 32*25);

      // Set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<2> cell;
      std::ifstream in;
      openInputFile("in/Rectangular", in);
      in >> cell;

      // set up wavelist object
      WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<2> const & ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<2> const & ksq_h2 = wavelist.kSq();

      // Compute kSq on host and compare
      IntVec<2> temp, vec;
      double ksq;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh2.dimensions(), cell);
         ksq = cell.ksq(vec);

         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(ksq - ksq_h2[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq3D()
   {
      printMethod(TEST_FUNC);

      // Set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<3> cell;
      std::ifstream in;
      openInputFile("in/Hexagonal", in);
      in >> cell;

      // Set up wavelist object
      WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq on device two different ways, transfer to host
      wavelist.computeMinimumImages(); // calculates kSq
      RField<3> const & ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<3> const & ksq_h2 = wavelist.kSq();

      // Compute kSq on host and compare
      IntVec<3> temp, vec;
      double ksq;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell);
         ksq = cell.ksq(vec);

         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(ksq - ksq_h2[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq1D()
   {
      printMethod(TEST_FUNC);

      // Set up wavelist object
      WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute dKSq on device, transfer to host
      wavelist.computedKSq();
      RField<1> dksq_h = wavelist.dKSq(0);

      // Compute dKSq on host and compare
      IntVec<1> temp, vec;
      double dksq;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         dksq = cell1.dksq(vec, 0);
         //if (mesh1.dimension(0) - temp[0] > mesh1.dimension(0)/2 + 1) {
         //   dksq *= 2;
         //}
         TEST_ASSERT(abs(dksq - dksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq2D()
   {
      printMethod(TEST_FUNC);

      WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);
      wavelist.computedKSq();

      // compute dKSq on host and compare
      IntVec<2> temp, vec;
      double dksq;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (int n = 0; n < cell2.nParameter() ; ++n) {
         RField<2> const & dksq_h = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            temp = iter.position();
            vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
            dksq = cell2.dksq(vec, n);
            //if (temp[1] != 0) {
               //if (mesh2.dimension(1) - temp[1] > mesh2.dimension(1)/2 + 1) {
               //   dksq *= 2;
               //}
            //}
            TEST_ASSERT(abs(dksq - dksq_h[iter.rank()]) < tolerance_);
         }
      }
   }

   void testComputedKSq3D()
   {
      printMethod(TEST_FUNC);

      WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq();

      IntVec<3> temp, vec;
      double dksq;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         RField<3> const & dksq_h = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            temp = iter.position();
            vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
            dksq = cell3.dksq(vec, n);
            //if (temp[2] != 0) {
            //   if (mesh3.dimension(2) - temp[2] > mesh3.dimension(2)/2 + 1) {
            //      dksq *= 2;
            //   }
            //}
            TEST_ASSERT(abs(dksq - dksq_h[iter.rank()]) < tolerance_);
         }
      }
   }

};

TEST_BEGIN(CpuWaveListTest)
TEST_ADD(CpuWaveListTest, testAllocate)
TEST_ADD(CpuWaveListTest, testComputeMinimumImages1D)
TEST_ADD(CpuWaveListTest, testComputeMinimumImages2D)
TEST_ADD(CpuWaveListTest, testComputeMinimumImages3D)
TEST_ADD(CpuWaveListTest, testComputeKSq1D)
TEST_ADD(CpuWaveListTest, testComputeKSq2D)
TEST_ADD(CpuWaveListTest, testComputeKSq3D)
TEST_ADD(CpuWaveListTest, testComputedKSq1D)
TEST_ADD(CpuWaveListTest, testComputedKSq2D)
TEST_ADD(CpuWaveListTest, testComputedKSq3D)
TEST_END(CpuWaveListTest)

#endif
