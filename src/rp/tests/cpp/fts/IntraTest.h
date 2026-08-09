#ifndef RPC_INTRA_TEST_H
#define RPC_INTRA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/compressor/IntraCorrelation.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldDft.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RFieldComparison.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/PolymerModel.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/random/Random.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class IntraTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   void testIntra(double A, double f, std::string paramFilename,
                  std::string inFieldFilename,
                  char const * outfilename)
   {
      openLogFile(outfilename);
      Rp::System<1,CPT> system;
      initSystem(system, paramFilename);
      system.w().readRGrid(inFieldFilename);

      int meshSize = system.domain().mesh().size();
      IntVec<1> const & dimensions = system.domain().mesh().dimensions();
      int nMonomer = system.mixture().nMonomer();
      double vMonomer = system.mixture().vMonomer();

      // Cos pressure field perturbation per chain: A * cos(2pi * f* i/meshSize)
      RField<1,CPT> cosF;
      RFieldDft<1,CPT> cosFK;
      cosF.allocate(dimensions);
      cosFK.allocate(dimensions);
      PolymerSpecies<double> const & polymer = system.mixture().polymerSpecies(0);
      for (int k = 0; k < meshSize; k++){
         cosF[k] = A * std::cos(2 * M_PI * k * f / meshSize);

         // Apply to each monomer
         if (PolymerModel::isBead()) {
            cosF[k] /= polymer.nBead();
         } else {
            cosF[k] /= polymer.length();
         }
      }

      // Convert to Fourier Space
      system.domain().fft().forwardTransform(cosF, cosFK);

      DArray< RField<1,CPT> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
         w2[i] = system.w().rgrid(i);
      }

      // Loop over monomer types add pressure perturbation
      for (int i = 0; i < nMonomer; ++i) {
         for (int k = 0; k < meshSize; ++k) {
            w2[i][k] += cosF[k];
         }
      }
      system.w().setRGrid(w2);

      system.compute();

      // Incompressibility error
      RField<1,CPT> error;
      error.allocate(dimensions);
      for (int k = 0; k <  meshSize; k++){
         error[k] = -1.0;
         for (int i = 0; i <nMonomer; i++){
            error[k] +=  system.c().rgrid(i)[k];
         }
      }

      // Intra analytical
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CPT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intra_(system);
      intra_.computeOmegaTotal(intraCorrelationK);

      // Compute analytical dphi using Intra
      RField<1,CPT> analyticalError;
      RFieldDft<1,CPT> analyticalErrorK;
      analyticalError.allocate(dimensions);
      analyticalErrorK.allocate(dimensions);
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDimensions);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         analyticalErrorK[iter.rank()][0] = -cosFK[iter.rank()][0] *  vMonomer * intraCorrelationK[iter.rank()];
         analyticalErrorK[iter.rank()][1] = -cosFK[iter.rank()][1] *  vMonomer * intraCorrelationK[iter.rank()];
      }

      system.domain().fft().inverseTransformUnsafe(analyticalErrorK, analyticalError);

      RFieldComparison<1,CPT> comparison;
      comparison.compare(error, analyticalError);
      TEST_ASSERT(comparison.maxDiff() < A* 1e-2);

   }

   void testIntraDiblockThread()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_diblock_thread",
               "in/w_diblock_homogenous.rf",
               "out/testIntraDiblockThread.log");

   }


   void testIntraTriblockThread()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_triblock_thread",
               "in/w_triblock_homogenous.rf",
               "out/testIntraTriblockThread.log");

   }

   void testIntraDiblockBead()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_diblock_bead",
               "in/w_diblock_homogenous.rf",
               "out/testIntraDiblockThread.log");

   }

   void testIntraTriblockBead()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_triblock_bead",
               "in/w_triblock_homogenous.rf",
               "out/testIntraTriblockBead.log");

   }

   void testIntraHomoThread()
   {
      // Compare the intracorrelation function of homopolymer and conformational diblock
      printMethod(TEST_FUNC);

      openLogFile("out/testIntraHomoThread.log");
      Rp::System<1,CPT> system;
      initSystem(system, "in/param_system_1D_diblock_thread");
      system.w().readRGrid("in/w_diblock_homogenous.rf");

      IntVec<1> const & dimensions = system.domain().mesh().dimensions();

      // The intracorrelation function of conformational diblock
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CPT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intra(system);
      intra.computeOmegaTotal(intraCorrelationK);

      // The intracorrelation function of conformational homo
      Rp::System<1,CPT> systemHomo;
      initSystem(systemHomo, "in/param_system_1D_homo_thread");
      systemHomo.w().readRGrid("in/w_homo_homogenous.rf");
      RField<1,CPT> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intraHomo(systemHomo);
      intraHomo.computeOmegaTotal(intraCorrelationKHomo);

      RFieldComparison<1,CPT> comparison;
      comparison.compare(intraCorrelationK, intraCorrelationKHomo);
      TEST_ASSERT(comparison.maxDiff() < 1e-5);

   }

   void testIntraHomoBead()
   {
      // Compare intracorrelation of homopolymer and conformational diblock
      printMethod(TEST_FUNC);

      openLogFile("out/testIntraHomoBead.log");
      Rp::System<1,CPT> system;
      initSystem(system, "in/param_system_1D_diblock_bead");
      system.w().readRGrid("in/w_diblock_homogenous.rf");

      IntVec<1> const & dimensions = system.domain().mesh().dimensions();

      // The intracorrelation function of conformational diblock
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CPT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intra(system);
      intra.computeOmegaTotal(intraCorrelationK);

      // The intracorrelation function of conformational homo
      Rp::System<1,CPT> systemHomo;
      initSystem(systemHomo, "in/param_system_1D_homo_bead");
      systemHomo.w().readRGrid("in/w_homo_homogenous.rf");
      RField<1,CPT> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intraHomo(systemHomo);
      intraHomo.computeOmegaTotal(intraCorrelationKHomo);

      RFieldComparison<1,CPT> comparison;
      comparison.compare(intraCorrelationK, intraCorrelationKHomo);
      TEST_ASSERT(comparison.maxDiff() < 1e-5);

   }

};

TEST_BEGIN(IntraTest)
TEST_ADD(IntraTest, testIntraDiblockThread)
TEST_ADD(IntraTest, testIntraTriblockThread)
TEST_ADD(IntraTest, testIntraDiblockBead)
TEST_ADD(IntraTest, testIntraTriblockBead)
TEST_ADD(IntraTest, testIntraHomoThread)
TEST_ADD(IntraTest, testIntraHomoBead)
TEST_END(IntraTest)

#endif
