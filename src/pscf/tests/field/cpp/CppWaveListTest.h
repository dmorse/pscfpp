#ifndef PSCF_CPP_WAVE_LIST_TEST_H
#define PSCF_CPP_WAVE_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/WaveList.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;

class CppWaveListTest : public UnitTest
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

      WaveList<1,CPT> wavelist1;
      wavelist1.allocate(mesh1, cell1);
      TEST_ASSERT(wavelist1.isAllocated());
      TEST_ASSERT(!wavelist1.hasMinImages());

      WaveList<2,CPT> wavelist2;
      wavelist2.allocate(mesh2, cell2);
      TEST_ASSERT(wavelist2.isAllocated());
      TEST_ASSERT(!wavelist2.hasMinImages());

      WaveList<3,CPT> wavelist3;
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
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<1> > const & minImages = wavelist.minImages();
      RField<1,CPT> const & ksq = wavelist.kSq();
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
         TEST_ASSERT((flag == FFT<1,CPT>::hasImplicitInverse(temp, meshDims1)));
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
      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell2);

      // compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<2> > const & minImages = wavelist.minImages();
      RField<2,CPT> ksq = wavelist.kSq();
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
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell3);

      // Compute minimum images (and ksq)
      wavelist.computeMinimumImages(); 
      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      RField<3,CPT> const & ksq = wavelist.kSq();

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
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<1,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq 
      RField<1,CPT> const & ksq2 = wavelist.kSq();

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
      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<2,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<2,CPT> const & ksq2 = wavelist.kSq();

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
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq two different ways
      wavelist.computeMinimumImages(); // calculates kSq
      RField<3,CPT> const & ksq = wavelist.kSq();
      wavelist.clearUnitCellData(); // resets kSq but not min images
      wavelist.computeKSq(); // recalculates kSq using a different kernel
      RField<3,CPT> const & ksq2 = wavelist.kSq();

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
      WaveList<1,CPT> wavelist;
      wavelist.allocate(mesh1, cell1);

      // Compute dKSq
      wavelist.computedKSq();
      RField<1,CPT> dksq = wavelist.dKSq(0);

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
            if (mesh1.dimension(0) - pos[0] > mesh1.dimension(0)/2) {
               val *= 2.0;
            }
         }
         TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
      }
   }

   void testComputedKSq2D()
   {
      printMethod(TEST_FUNC);

      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell2);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<2> pos, vec;
      double val;
      MeshIterator<2> iter;
      iter.setDimensions(kMeshDims2);
      for (int n = 0; n < cell2.nParameter() ; ++n) {
         RField<2,CPT> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh2.dimensions(), cell2);
            val = cell2.dksq(vec, n);
            if (pos[1] != 0) {
               if (mesh2.dimension(1) - pos[1] > mesh2.dimension(1)/2){
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

      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq();

      // compute dKSq locally and compare
      IntVec<3> pos, vec;
      double val;
      MeshIterator<3> iter;
      iter.setDimensions(kMeshDims3);
      for (int n = 0; n < cell3.nParameter() ; ++n) {
         RField<3,CPT> const & dksq = wavelist.dKSq(n);
         for (iter.begin(); !iter.atEnd(); ++iter) {
            pos = iter.position();
            vec = shiftToMinimum(pos, mesh3.dimensions(), cell3);
            val = cell3.dksq(vec, n);
            if (pos[2] != 0) {
               if (mesh3.dimension(2) - pos[2] > mesh3.dimension(2)/2){
                  val *= 2.0;
               }
            }
            TEST_ASSERT(abs(val - dksq[iter.rank()]) < tolerance_);
         }
      }
   }

   void testSortWaves3D()
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
      WaveList<3,CPT> wavelist;
      wavelist.allocate(mesh3, cell);

      // Compute kSq 
      wavelist.computeMinimumImages(); 
      wavelist.computeKSq(); 
      RField<3,CPT> const & ksq = wavelist.kSq();
      int kSize = wavelist.kSize();
      TEST_ASSERT(kSize > 0);

      wavelist.sortWaves();
      int nBunch = wavelist.nBunch(); 
     
      // Test sortedIds 
      DArray<int> const & sortedIds = wavelist.sortedIds();
      for (int i = 1; i < kSize; ++i) {
         TEST_ASSERT(ksq[sortedIds[i]] >= ksq[sortedIds[i-1]]);
      }

      // Test sortedBunches
      double epsilon = 1.0E-8;
      double value = 0.0;
      double beginValue = 0.0;
      int begin, end, size, ib, iw;
      GArray< Pair<int> > const & sortedBunches = wavelist.sortedBunches();
      for (ib = 0; ib < nBunch; ++ib) {
         begin = sortedBunches[ib][0];
         end = sortedBunches[ib][1];
         size = end - begin;
         value = ksq[sortedIds[begin]];
         if (ib > 0) {
            TEST_ASSERT(std::abs(value - beginValue) > epsilon);
         }
         beginValue = value;
         if (size > 1) {
            for (iw = begin + 1; iw < end; ++iw) {
               value = ksq[sortedIds[iw]];
               TEST_ASSERT(std::abs(value - beginValue) < 2.0*epsilon);
            }
         }
      }

      // Test bunchIds
      TEST_ASSERT(nBunch > 0);
      DArray<int> const & bunchIds = wavelist.bunchIds();
      TEST_ASSERT(bunchIds.capacity() == kSize);
      DArray<double> oldVal(nBunch);
      DArray<bool> found(nBunch);
      double newVal;
      for (ib = 0; ib < nBunch; ++ib) {
         found[ib] = false;
      }
      for (iw = 0; iw < kSize; ++iw) {
         newVal = ksq[iw];
         ib = bunchIds[iw];
         if (!found[ib]) {
            oldVal[ib] = newVal;
            found[ib] = true;
         } else {
            TEST_ASSERT(std::abs(newVal - oldVal[ib]) < 2.0E-8);
         }
      }
   }

   void testComplex()
   {
      printMethod(TEST_FUNC);

      WaveList<3,CPT> wavelist(false);
      wavelist.allocate(mesh3, cell3);
      wavelist.computedKSq(); // computes min images, ksq, and dksq

      DArray< IntVec<3> > const & minImages = wavelist.minImages();
      RField<3,CPT> const & ksq = wavelist.kSq();
      DArray< RField<3,CPT> > const & dksq = wavelist.dKSq();

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

   /*
   * Sums S[n] over the k-grid of w dKSq_n, with w = exp(-kSq/kSqMax),
   * for a cell and a mesh with the symmetry of a hexagonal or cubic
   * lattice. The tensor sum of w k k is then isotropic, so that S[n]
   * is proportional to dLogV[n], the derivative of the logarithm of
   * the cell volume. Returns the largest deviation from this relation,
   * relative to the largest |S[n]|.
   */
   template <int D>
   double dKSqAnisotropy(UnitCell<D> const & cell, IntVec<D> const & d,
                         double const * dLogV)
   {
      Mesh<D> mesh(d);
      WaveList<D,CPT> wavelist;
      wavelist.allocate(mesh, cell);
      wavelist.computedKSq();
      RField<D,CPT> const & kSq = wavelist.kSq();
      double kSqMax = 0.0;
      for (int i = 0; i < kSq.capacity(); ++i) {
         kSqMax = std::max(kSqMax, kSq[i]);
      }
      double S[6];
      double scale = 0.0;
      for (int n = 0; n < cell.nParameter(); ++n) {
         S[n] = 0.0;
         for (int i = 0; i < kSq.capacity(); ++i) {
            S[n] += std::exp(-kSq[i]/kSqMax)*wavelist.dKSq(n)[i];
         }
         scale = std::max(scale, std::abs(S[n]));
      }
      double deviation = 0.0;
      for (int n = 1; n < cell.nParameter(); ++n) {
         double r = S[n] - S[0]*dLogV[n]/dLogV[0];
         deviation = std::max(deviation, std::abs(r)/scale);
      }
      return deviation;
   }

   /*
   * Cells with a variable angle at geometries of higher symmetry, in
   * which the minimum images of some wavevectors are tied with images
   * whose |k|^2 has other derivatives: dKSq must treat these images
   * alike, so that the sums of dKSqAnisotropy are isotropic. Elsewhere,
   * dKSq must be the derivative for the minimum image, bit for bit.
   */
   void testComputedKSqTied()
   {
      printMethod(TEST_FUNC);

      // Rhombic cell at 120 degrees (hexagonal lattice), area a^2 sin(g)
      FSArray<double, 6> parameters;
      double a = 1.7;
      double g = 2.0*Constants::Pi/3.0;
      parameters.append(a);
      parameters.append(g);
      UnitCell<2> rhombic;
      rhombic.set(UnitCell<2>::Rhombic, parameters);
      double dLogV[2];
      dLogV[0] = 2.0/a;
      dLogV[1] = std::cos(g)/std::sin(g);
      TEST_ASSERT(dKSqAnisotropy<2>(rhombic, IntVec<2>(12), dLogV) < 1E-12);
      TEST_ASSERT(dKSqAnisotropy<2>(rhombic, IntVec<2>(13), dLogV) < 1E-12);

      // Rhombohedral cells at 60 degrees (fcc) and 109.47 degrees (bcc),
      // volume a^3 (1 - 3c^2 + 2c^3)^(1/2), with c = cos(alpha)
      a = 2.0;
      double alpha[2] = {Constants::Pi/3.0, std::acos(-1.0/3.0)};
      for (int k = 0; k < 2; ++k) {
         parameters[0] = a;
         parameters[1] = alpha[k];
         UnitCell<3> rhombohedral;
         rhombohedral.set(UnitCell<3>::Rhombohedral, parameters);
         double c = std::cos(alpha[k]);
         dLogV[0] = 3.0/a;
         dLogV[1] = 3.0*c*(1.0 - c)*std::sin(alpha[k])
                    /(1.0 - 3.0*c*c + 2.0*c*c*c);
         TEST_ASSERT(dKSqAnisotropy<3>(rhombohedral, IntVec<3>(8), dLogV)
                     < 1E-12);
      }

      // Generic oblique cell: dKSq of the minimum image, bit for bit
      WaveList<2,CPT> wavelist;
      wavelist.allocate(mesh2, cell2);
      wavelist.computedKSq();
      for (int n = 0; n < cell2.nParameter(); ++n) {
         for (int i = 0; i < kSize2; ++i) {
            double dksq = cell2.dksq(wavelist.minImages()[i], n);
            if (wavelist.implicitInverse()[i]) dksq *= 2.0;
            TEST_ASSERT(wavelist.dKSq(n)[i] == dksq);
         }
      }
   }

};

TEST_BEGIN(CppWaveListTest)
TEST_ADD(CppWaveListTest, testAllocate)
TEST_ADD(CppWaveListTest, testComputeMinimumImages1D)
TEST_ADD(CppWaveListTest, testComputeMinimumImages2D)
TEST_ADD(CppWaveListTest, testComputeMinimumImages3D)
TEST_ADD(CppWaveListTest, testComputeKSq1D)
TEST_ADD(CppWaveListTest, testComputeKSq2D)
TEST_ADD(CppWaveListTest, testComputeKSq3D)
TEST_ADD(CppWaveListTest, testComputedKSq1D)
TEST_ADD(CppWaveListTest, testComputedKSq2D)
TEST_ADD(CppWaveListTest, testComputedKSq3D)
TEST_ADD(CppWaveListTest, testSortWaves3D)
TEST_ADD(CppWaveListTest, testComplex)
TEST_ADD(CppWaveListTest, testComputedKSqTied)
TEST_END(CppWaveListTest)

#endif
