#ifndef PRDC_CPU_WAVE_LIST_TEST_H
#define PRDC_CPU_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/WaveList.h>
#include <prdc/cpu/RField.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

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

   IntVec<1> meshDims1;
   IntVec<2> meshDims2;
   IntVec<3> meshDims3;

   IntVec<1> kMeshDims1;
   IntVec<2> kMeshDims2;
   IntVec<3> kMeshDims3;

   int kSize1, kSize2, kSize3;

public:

   void setUp()
   {
      // set up 1D mesh
      // IntVec<1> meshDims1; 
      meshDims1[0] = 32; 
      mesh1.setDimensions(meshDims1);
      kMeshDims1[0] = (meshDims1[0] / 2) + 1;
      kSize1 = kMeshDims1[0];

      // Set up 2D meshes
      //IntVec<2> meshDims2; 
      meshDims2[0] = 32; 
      meshDims2[1] = 48; 
      mesh2.setDimensions(meshDims2);
      kMeshDims2 = meshDims2;
      kMeshDims2[1] = meshDims2[1] / 2 + 1;
      kSize2 = kMeshDims2[0] * kMeshDims2[1];

      // Set up 3D meshes
      //IntVec<3> meshDims3; 
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

      Cpu::WaveList<1> wavelist1;
      wavelist1.allocate(mesh1, cell1);
      TEST_ASSERT(wavelist1.isAllocated());
      TEST_ASSERT(!wavelist1.hasMinImages());

      Cpu::WaveList<2> wavelist2;
      wavelist2.allocate(mesh2, cell2);
      TEST_ASSERT(wavelist2.isAllocated());
      TEST_ASSERT(!wavelist2.hasMinImages());

      Cpu::WaveList<3> wavelist3;
      wavelist3.allocate(mesh3, cell3);
      TEST_ASSERT(wavelist3.isAllocated());
      TEST_ASSERT(!wavelist3.hasMinImages());
   }

   void testComputeMinimumImages1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] == 17);
      TEST_ASSERT(kSize1 == 17);

      // set up wavelist object
      Cpu::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<1> > const & minImages = wavelist.minImages();
      Cpu::RField<1> const & ksq = wavelist.kSq();
      DArray<bool> const & implicit = wavelist.implicitInverse();

      // Compute minimum images (and ksq) locally and compare
      IntVec<1> temp, vec;
      MeshIterator<1> iter;
      double val;
      int rank;
      bool flag;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         flag = implicit[rank];
         val = cell1.ksq(vec);
         
         TEST_ASSERT(vec == minImages[rank]);
         TEST_ASSERT(abs(val - ksq[rank]) < tolerance_);
         TEST_ASSERT(flag == Cpu::FFT<1>::hasImplicitInverse(temp, meshDims1));
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
      Cpu::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<2> > const & minImages = wavelist.minImages();
      Cpu::RField<2> ksq = wavelist.kSq();
      TEST_ASSERT(minImages.capacity() == kSize2);
      TEST_ASSERT(ksq.capacity() == kSize2);

      // compute minimum images (and ksq) locally and compare
      IntVec<2> temp, vec;
      MeshIterator<2> iter;
      double val;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
         val = cell2.ksq(vec);
         
         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
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
      Cpu::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      Cpu::RField<3> const & ksq = wavelist.kSq();

      // Compute minimum images (and ksq) locally and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double val;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         val = cell3.ksq(vec);

         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq1D()
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(mesh1.dimension(0) == 32);
      TEST_ASSERT(kMeshDims1[0] == 17);
      TEST_ASSERT(kSize1 == 17);

      // set up wavelist object
      Cpu::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      Cpu::RField<1> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq 
      Cpu::RField<1> const & ksq2 = wavelist.kSq();

      // Compute kSq locally and compare
      IntVec<1> temp, vec;
      double val;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         val = cell1.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
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

      // Set up wavelist object
      Cpu::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      Cpu::RField<2> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      Cpu::RField<2> const & ksq2 = wavelist.kSq();

      // Compute kSq in wavelist and test
      IntVec<2> pos, vec;
      double val;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         pos = iter.position();
         vec = shiftToMinimum(pos, mesh2.dimensions(), cell);
         val = cell.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
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
      Cpu::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      Cpu::RField<3> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      Cpu::RField<3> const & ksq2 = wavelist.kSq();

      // Compute kSq locally and compare
      IntVec<3> temp, vec;
      double val;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell);
         val = cell.ksq(vec);

         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);
         TEST_ASSERT(abs(val - ksq2[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq1D()
   {
      printMethod(TEST_FUNC);

      // Set up wavelist object
      Cpu::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute dKSq
      wavelist.computedKSq();
      Cpu::RField<1> dksq = wavelist.dKSq(0);

      // Compute dKSq locally and compare
      IntVec<1> pos, vec;
      double val;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         pos = iter.position();
         vec = shiftToMinimum(pos, mesh1.dimensions(), cell1);
         val = cell1.dksq(vec, 0);
         if (pos[0] != 0) {
            if (mesh1.dimension(0) - pos[0] > mesh1.dimension(0)/2 + 1) {
               val *= 2.0;
            }
         }
         TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq2D()
   {
      printMethod(TEST_FUNC);

      Cpu::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<2> pos, vec;
      double val;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (int n = 0; n < cell2.nParameter() ; ++n) {
         Cpu::RField<2> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh2.dimensions(), cell2);
            val = cell2.dksq(vec, n);
            if (pos[1] != 0) {
               if (mesh2.dimension(1)-pos[1] > mesh2.dimension(1)/2 + 1){
                  val *= 2.0;
               }
            }
            TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
         }
      }
   }

   void testComputedKSq3D()
   {
      printMethod(TEST_FUNC);

      Cpu::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<3> pos, vec;
      double val;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         Cpu::RField<3> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh3.dimensions(), cell3);
            val = cell3.dksq(vec, n);
            if (pos[2] != 0) {
               if (mesh3.dimension(2)-pos[2] > mesh3.dimension(2)/2 + 1){
                  val *= 2.0;
               }
            }
            TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
         }
      }
   }

   void testComplex()
   {
      printMethod(TEST_FUNC);

      Cpu::WaveList<3> wavelist(false);
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq(); // computes min images, ksq, and dksq

      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      Cpu::RField<3> const & ksq = wavelist.kSq();
      DArray< Cpu::RField<3> > const & dksq = wavelist.dKSq();

      // Check that array sizes are correct
      TEST_ASSERT(minImages.capacity() == mesh3.size());
      TEST_ASSERT(ksq.capacity() == mesh3.size());
      TEST_ASSERT(dksq.capacity() == cell3.nParameter());
      for (int i = 0; i < cell3.nParameter(); i++) {
         TEST_ASSERT(dksq[i].capacity() == mesh3.size());
      }

      // Compute minimum images, ksq, and dksq locally and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double val;
      iter.setDimensions(mesh3.dimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         val = cell3.ksq(vec);

         TEST_ASSERT(vec == minImages[iter.rank()]);
         TEST_ASSERT(abs(val - ksq[iter.rank()]) < tolerance_);

         for (int i = 0; i < cell3.nParameter(); i++) {
            val = cell3.dksq(vec, i);
            TEST_ASSERT(abs(val - dksq[i][iter.rank()]) < tolerance_);
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
TEST_ADD(CpuWaveListTest, testComplex)
TEST_END(CpuWaveListTest)

#endif
