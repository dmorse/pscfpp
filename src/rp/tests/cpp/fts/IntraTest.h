#ifndef RPC_INTRA_TEST_H
#define RPC_INTRA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/compressor/IntraCorrelation.h>
#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/PolymerModel.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/random/Random.h>

#include <fstream>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

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

   void testIntraThreadResponse()
   {
      // Compare the correlation of the thread model with the response
      // of the MDE solution, at all wavenumbers of a mesh on which
      // kSq*b*b*ds/6 is much larger than 1 at the largest ones
      printMethod(TEST_FUNC);
      openLogFile("out/testIntraThreadResponse.log");

      // Triblock with two segment lengths, and a solvent. On 16 grid
      // points over a length 0.5, kSq*b*b*ds/6 = 16.8 (b = 1) and 37.9
      // (b = 1.5) at the largest wavenumber
      Rp::System<1,CPT> system;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.5 ]"
         " nPolymer 1 nSolvent 1"
         " Polymer{ type linear nBlock 3 blocks[ 0 0.25 1 0.5 0 0.25 ]"
         " phi 0.8 }"
         " Solvent{ monomerId 1 size 0.5 phi 0.2 }"
         " vMonomer 0.01 ds 0.01 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar } }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(0.5);
      system.setUnitCell(parameters);
      const int meshSize = 16;
      const int nMonomer = 2;
      const double vMonomer = 0.01;

      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = meshSize/2 + 1;
      RField<1,CPT> omega;
      omega.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CPT> intra(system);
      intra.computeOmegaTotal(omega);

      // If all w fields change from 0 to eps*cos(2 pi n x/L), the total
      // concentration changes by -eps*vMonomer*omega*cos(2 pi n x/L), to
      // first order in eps
      const double eps = 1.0E-4;
      DArray< RField<1,CPT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
      }
      double maxDev = 0.0;
      for (int n = 1; n <= meshSize/2; ++n) {
         for (int j = 0; j < meshSize; ++j) {
            for (int i = 0; i < nMonomer; ++i) {
               w[i][j] = eps*std::cos(2.0*M_PI*n*j/meshSize);
            }
         }
         system.w().setRGrid(w);
         system.compute();
         double dphi = 0.0;
         for (int j = 0; j < meshSize; ++j) {
            double phi = 0.0;
            for (int i = 0; i < nMonomer; ++i) {
               phi += system.c().rgrid(i)[j];
            }
            dphi += phi*std::cos(2.0*M_PI*n*j/meshSize);
         }
         dphi *= (2*n == meshSize) ? 1.0/meshSize : 2.0/meshSize;
         double response = -dphi/(eps*vMonomer);
         Log::file() << "n = " << n << "  response = " << response
                     << "  omega = " << omega[n] << std::endl;
         double dev = std::abs(response/omega[n] - 1.0);
         TEST_ASSERT(dev < 1.0E-5);
         if (dev > maxDev) maxDev = dev;
      }
      Log::file() << "Largest relative deviation " << maxDev << std::endl;
   }

};

TEST_BEGIN(IntraTest)
TEST_ADD(IntraTest, testIntraDiblockThread)
TEST_ADD(IntraTest, testIntraTriblockThread)
TEST_ADD(IntraTest, testIntraDiblockBead)
TEST_ADD(IntraTest, testIntraTriblockBead)
TEST_ADD(IntraTest, testIntraHomoThread)
TEST_ADD(IntraTest, testIntraHomoBead)
TEST_ADD(IntraTest, testIntraThreadResponse)
TEST_END(IntraTest)

#endif
