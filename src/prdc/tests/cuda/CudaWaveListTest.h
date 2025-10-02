#ifndef PRDC_CUDA_WAVE_LIST_TEST_H
#define PRDC_CUDA_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/WaveList.h>
#include <prdc/cuda/types.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/cuda/HostDArray.h> 
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaWaveListTest : public UnitTest
{

//using namespace Pscf::Prdc::Cuda;

private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   constexpr static Cuda::cudaReal tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static Cuda::cudaReal tolerance_ = 1E-10;
   #endif
   #endif

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
      // set up 1D, 2D, and 3D meshes
      IntVec<1> meshDims1; 
      meshDims1[0] = 32; 
      mesh1.setDimensions(meshDims1);
      kMeshDims1[0] = meshDims1[0] / 2 + 1;
      kSize1 = kMeshDims1[0];

      IntVec<2> meshDims2; 
      meshDims2[0] = 32; 
      meshDims2[1] = 48; 
      mesh2.setDimensions(meshDims2);
      kMeshDims2 = meshDims2;
      kMeshDims2[1] = meshDims2[1] / 2 + 1;
      kSize2 = kMeshDims2[0] * kMeshDims2[1];

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

      using namespace Pscf::Prdc::Cuda;

      Cuda::WaveList<1> wavelist1;
      wavelist1.allocate(mesh1, cell1);
      TEST_ASSERT(wavelist1.isAllocated());
      TEST_ASSERT(!wavelist1.hasMinImages());

      Cuda::WaveList<2> wavelist2;
      wavelist2.allocate(mesh2, cell2);
      TEST_ASSERT(wavelist2.isAllocated());
      TEST_ASSERT(!wavelist2.hasMinImages());

      Cuda::WaveList<3> wavelist3;
      wavelist3.allocate(mesh3, cell3);
      TEST_ASSERT(wavelist3.isAllocated());
      TEST_ASSERT(!wavelist3.hasMinImages());
   }

   void testComputeMinimumImages1D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // compute minimum images (and ksq) on device, transfer to host
      HostDArray<int> minImages_h;
      HostDArray<Cuda::cudaReal> ksq_h;
      wavelist.computeMinimumImages(); 
      minImages_h = wavelist.minImages_d();
      ksq_h = wavelist.kSq();

      // compute minimum images (and ksq) on host and compare
      IntVec<1> temp, vec;
      MeshIterator<1> iter;
      Cuda::cudaReal ksq;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         ksq = cell1.ksq(vec);
         
         TEST_ASSERT(vec[0] == minImages_h[iter.rank()]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputeMinimumImages2D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute minimum images (and ksq) on device, transfer to host
      HostDArray<int> minImages_h;
      HostDArray<Cuda::cudaReal> ksq_h;
      wavelist.computeMinimumImages(); 
      minImages_h = wavelist.minImages_d();
      ksq_h = wavelist.kSq();

      // compute minimum images (and ksq) on host and compare
      IntVec<2> temp, vec;
      MeshIterator<2> iter;
      Cuda::cudaReal ksq;
      iter.setDimensions(kMeshDims2);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
         ksq = cell2.ksq(vec);
         
         TEST_ASSERT(vec[0] == minImages_h[iter.rank()]);
         TEST_ASSERT(vec[1] == minImages_h[iter.rank() + kSize2]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputeMinimumImages3D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);

      // compute minimum images (and ksq) on device, transfer to host
      HostDArray<int> minImages_h;
      HostDArray<Cuda::cudaReal> ksq_h;
      wavelist.computeMinimumImages(); 
      minImages_h = wavelist.minImages_d();
      ksq_h = wavelist.kSq();

      // compute minimum images (and ksq) on host and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      Cuda::cudaReal ksq;
      iter.setDimensions(kMeshDims3);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         ksq = cell3.ksq(vec);

         TEST_ASSERT(vec[0] == minImages_h[iter.rank()]);
         TEST_ASSERT(vec[1] == minImages_h[iter.rank() + kSize3]);
         TEST_ASSERT(vec[2] == minImages_h[iter.rank() + kSize3 + kSize3]);
         TEST_ASSERT(abs(ksq - ksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputeKSq1D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // compute kSq on device two different ways, transfer to host
      HostDArray<Cuda::cudaReal> ksq_h, ksq_h2;
      wavelist.computeMinimumImages(); // calculates kSq
      ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      ksq_h2 = wavelist.kSq();

      // compute kSq on host and compare
      IntVec<1> temp, vec;
      Cuda::cudaReal ksq;
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
      using namespace Pscf::Prdc::Cuda;

      // set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<2> cell;
      std::ifstream in;
      openInputFile("in/Rectangular", in);
      in >> cell;

      // set up wavelist object
      Cuda::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell);

      // compute kSq on device two different ways, transfer to host
      HostDArray<Cuda::cudaReal> ksq_h, ksq_h2;
      wavelist.computeMinimumImages(); // calculates kSq
      ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      ksq_h2 = wavelist.kSq();

      // compute kSq on host and compare
      IntVec<2> temp, vec;
      Cuda::cudaReal ksq;
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
      using namespace Pscf::Prdc::Cuda;

      // set up unit cell with no flexible angles
      // (if there are flexible angles, computeKSq never gets to run,
      //  because computeMinimumImages is always called instead.)
      UnitCell<3> cell;
      std::ifstream in;
      openInputFile("in/Hexagonal", in);
      in >> cell;

      // set up wavelist object
      Cuda::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell);

      // compute kSq on device two different ways, transfer to host
      HostDArray<Cuda::cudaReal> ksq_h, ksq_h2;
      wavelist.computeMinimumImages(); // calculates kSq
      ksq_h = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      ksq_h2 = wavelist.kSq();

      // compute kSq on host and compare
      IntVec<3> temp, vec;
      Cuda::cudaReal ksq;
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
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<1> wavelist;
      wavelist.allocate(mesh1, cell1);

      // compute dKSq on device, transfer to host
      wavelist.computedKSq();
      HostDArray<Cuda::cudaReal> dksq_h;
      dksq_h = wavelist.dKSq(0);

      // compute dKSq on host and compare
      IntVec<1> temp, vec;
      Cuda::cudaReal dksq;
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDims1);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh1.dimensions(), cell1);
         dksq = cell1.dksq(vec, 0);
         if (mesh1.dimension(0) - temp[0] > mesh1.dimension(0)/2 + 1) {
            dksq *= 2;
         }
         TEST_ASSERT(abs(dksq - dksq_h[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq2D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<2> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute dKSq on device, transfer to host
      wavelist.computedKSq();
      DArray< HostDArray<Cuda::cudaReal> > dksq_h;
      dksq_h.allocate(cell2.nParameter());
      for (int n = 0; n < cell2.nParameter() ; ++n) {
        dksq_h[n] = wavelist.dKSq(n);
      }

      // compute dKSq on host and compare
      IntVec<2> temp, vec;
      Cuda::cudaReal dksq;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (int n = 0; n < cell2.nParameter() ; ++n) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            temp = iter.position();
            vec = shiftToMinimum(temp, mesh2.dimensions(), cell2);
            dksq = cell2.dksq(vec, n);
            if (temp[1] != 0) {
               if (mesh2.dimension(1) - temp[1] > mesh2.dimension(1)/2 + 1) {
                  dksq *= 2;
               }
            }
            TEST_ASSERT(abs(dksq - dksq_h[n][iter.rank()]) < tolerance_);
         }
      }
   }

   void testComputedKSq3D()
   {
      printMethod(TEST_FUNC);
      using namespace Pscf::Prdc::Cuda;

      // set up wavelist object
      Cuda::WaveList<3> wavelist;
      wavelist.allocate(mesh3, cell3);

      // compute dKSq on device, transfer to host
      wavelist.computedKSq();
      DArray< HostDArray<Cuda::cudaReal> > dksq_h;
      dksq_h.allocate(cell3.nParameter());
      for (int n = 0; n < cell3.nParameter() ; ++n) {
        dksq_h[n] = wavelist.dKSq(n);
      }

      // compute dKSq on host and compare
      IntVec<3> temp, vec;
      Cuda::cudaReal dksq;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            temp = iter.position();
            vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
            dksq = cell3.dksq(vec, n);
            if (temp[2] != 0) {
               if (mesh3.dimension(2) - temp[2] > mesh3.dimension(2)/2 + 1) {
                  dksq *= 2;
               }
            }
            TEST_ASSERT(abs(dksq - dksq_h[n][iter.rank()]) < tolerance_);
         }
      }
   }

   void testComplex()
   {
      printMethod(TEST_FUNC);

      Cuda::WaveList<3> wavelist(false);
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq(); // computes min images, ksq, and dksq

      // Transfer results to host
      HostDArray<int> minImages_h;
      minImages_h = wavelist.minImages_d();

      HostDArray<Cuda::cudaReal> ksq_h;
      ksq_h = wavelist.kSq();

      DArray< HostDArray<Cuda::cudaReal> > dksq_h;
      dksq_h.allocate(cell3.nParameter());
      for (int n = 0; n < cell3.nParameter() ; ++n) {
        dksq_h[n] = wavelist.dKSq(n);
      }

      // Check that array sizes are correct
      int meshSize = mesh3.size();
      TEST_ASSERT(minImages_h.capacity() == meshSize * 3);
      TEST_ASSERT(ksq_h.capacity() == meshSize);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         TEST_ASSERT(dksq_h[n].capacity() == meshSize);
      }

      // Compute minimum images, ksq, and dksq on host and compare
      IntVec<3> temp, vec;
      MeshIterator<3> iter;
      double val;
      iter.setDimensions(mesh3.dimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         temp = iter.position();
         vec = shiftToMinimum(temp, mesh3.dimensions(), cell3);
         val = cell3.ksq(vec);

         TEST_ASSERT(vec[0] == minImages_h[iter.rank()]);
         TEST_ASSERT(vec[1] == minImages_h[iter.rank() + meshSize]);
         TEST_ASSERT(vec[2] == minImages_h[iter.rank() + meshSize + meshSize]);
         TEST_ASSERT(abs(val - ksq_h[iter.rank()]) < tolerance_);

         for (int n = 0; n < cell3.nParameter(); n++) {
            val = cell3.dksq(vec, n);
            TEST_ASSERT(abs(val - dksq_h[n][iter.rank()]) < tolerance_);
         }
      }
   }

};

TEST_BEGIN(CudaWaveListTest)
TEST_ADD(CudaWaveListTest, testAllocate)
TEST_ADD(CudaWaveListTest, testComputeMinimumImages1D)
TEST_ADD(CudaWaveListTest, testComputeMinimumImages2D)
TEST_ADD(CudaWaveListTest, testComputeMinimumImages3D)
TEST_ADD(CudaWaveListTest, testComputeKSq1D)
TEST_ADD(CudaWaveListTest, testComputeKSq2D)
TEST_ADD(CudaWaveListTest, testComputeKSq3D)
TEST_ADD(CudaWaveListTest, testComputedKSq1D)
TEST_ADD(CudaWaveListTest, testComputedKSq2D)
TEST_ADD(CudaWaveListTest, testComputedKSq3D)
TEST_ADD(CudaWaveListTest, testComplex)
TEST_END(CudaWaveListTest)

#endif
