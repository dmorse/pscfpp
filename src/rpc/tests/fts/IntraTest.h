#ifndef RPC_INTRA_TEST_H
#define RPC_INTRA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/intra/IntraCorrelation.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/PolymerModel.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/random/Random.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

class IntraTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(System<D>& system, std::string filename)
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
      System<1> system;
      initSystem(system, paramFilename);
      system.readWRGrid(inFieldFilename);

      int meshSize = system.domain().mesh().size();
      IntVec<1> const & dimensions = system.domain().mesh().dimensions();
      int nMonomer = system.mixture().nMonomer();
      double vMonomer = system.mixture().vMonomer();

      // Cos pressure field perturbation per chain: A * cos(2pi * f* i/meshSize)
      RField<1> cosF;
      RFieldDft<1> cosFK;
      cosF.allocate(dimensions);
      cosFK.allocate(dimensions);
      PolymerSpecies const & polymer = system.mixture().polymerSpecies(0);
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

      DArray< RField<1> > w2;
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
      system.setWRGrid(w2);

      system.compute();

      // Incompressibility error
      RField<1> error;
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
      RField<1> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      IntraCorrelation<1> intra_(system);
      intra_.computeIntraCorrelations(intraCorrelationK);

      // Compute analytical dphi using Intra
      RField<1> analyticalError;
      RFieldDft<1> analyticalErrorK;
      analyticalError.allocate(dimensions);
      analyticalErrorK.allocate(dimensions);
      MeshIterator<1> iter;
      iter.setDimensions(kMeshDimensions);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         analyticalErrorK[iter.rank()][0] = -cosFK[iter.rank()][0] *  vMonomer * intraCorrelationK[iter.rank()];
         analyticalErrorK[iter.rank()][1] = -cosFK[iter.rank()][1] *  vMonomer * intraCorrelationK[iter.rank()];
      }

      system.domain().fft().inverseTransformUnsafe(analyticalErrorK, analyticalError);

      RFieldComparison<1> comparison;
      comparison.compare(error, analyticalError);
      TEST_ASSERT(comparison.maxDiff() < A* 1e-2);

   }

   void testIntraDiblockThread()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_diblcok_thread",
               "in/w_diblock_homogenous.rf",
               "out/testIntraDiblockThread.log");

   }


   void testIntraTriblockThread()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_triblcok_thread",
               "in/w_triblock_homogenous.rf",
               "out/testIntraTriblockThread.log");

   }

   void testIntraDiblockBead()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_diblcok_bead",
               "in/w_diblock_homogenous.rf",
               "out/testIntraDiblockThread.log");

   }

   void testIntraTriblockBead()
   {
      printMethod(TEST_FUNC);
      double A = 1e-3;
      double f = 1.0;
      testIntra(A, f, "in/param_system_1D_triblcok_bead",
               "in/w_triblock_homogenous.rf",
               "out/testIntraTriblockBead.log");

   }

   void testIntraHomoThread()
   {
      // Compare the intracorrelation function of homopolymer and conformational diblock
      printMethod(TEST_FUNC);

      openLogFile("out/testIntraHomoThread.log");
      System<1> system;
      initSystem(system, "in/param_system_1D_diblcok_thread");
      system.readWRGrid("in/w_diblock_homogenous.rf");

      IntVec<1> const & dimensions = system.domain().mesh().dimensions();

      // The intracorrelation function of conformational diblock
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      IntraCorrelation<1> intra(system);
      intra.computeIntraCorrelations(intraCorrelationK);

      // The intracorrelation function of conformational homo
      System<1> systemHomo;
      initSystem(systemHomo, "in/param_system_1D_homo_thread");
      systemHomo.readWRGrid("in/w_homo_homogenous.rf");
      RField<1> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      IntraCorrelation<1> intraHomo(systemHomo);
      intraHomo.computeIntraCorrelations(intraCorrelationKHomo);

      RFieldComparison<1> comparison;
      comparison.compare(intraCorrelationK, intraCorrelationKHomo);
      TEST_ASSERT(comparison.maxDiff() < 1e-5);

   }

   void testIntraHomoBead()
   {
      // Compare intracorrelation of homopolymer and conformational diblock
      printMethod(TEST_FUNC);

      openLogFile("out/testIntraHomoBead.log");
      System<1> system;
      initSystem(system, "in/param_system_1D_diblcok_bead");
      system.readWRGrid("in/w_diblock_homogenous.rf");

      IntVec<1> const & dimensions = system.domain().mesh().dimensions();

      // The intracorrelation function of conformational diblock
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      IntraCorrelation<1> intra(system);
      intra.computeIntraCorrelations(intraCorrelationK);

      // The intracorrelation function of conformational homo
      System<1> systemHomo;
      initSystem(systemHomo, "in/param_system_1D_homo_bead");
      systemHomo.readWRGrid("in/w_homo_homogenous.rf");
      RField<1> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      IntraCorrelation<1> intraHomo(systemHomo);
      intraHomo.computeIntraCorrelations(intraCorrelationKHomo);

      RFieldComparison<1> comparison;
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
