#ifndef RPG_INTRA_TEST_H
#define RPG_INTRA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/fts/compressor/IntraCorrelation.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/field/CFields.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/RField.h>
#include <prdc/field/cuda/RFieldDft.h>
#include <prdc/field/cuda/RFieldComparison.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/PolymerModel.h>

#include <util/tests/LogFileUnitTest.h>

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
   void initSystem(Rp::System<D,CUT>& system, std::string filename)
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
      Rp::System<1,CUT> system; 
      initSystem(system, paramFilename);
      system.w().readRGrid(inFieldFilename);
      
      int meshSize = system.domain().mesh().size();
      IntVec<1> const & dimensions = system.domain().mesh().dimensions();
      int nMonomer = system.mixture().nMonomer();
      double vMonomer = system.mixture().vMonomer();
      
      // Cos pressure field perturbation per chain: A * cos(2pi * f* i/meshSize)
      RField<1,CUT> cosF;
      RFieldDft<1,CUT> cosFK;
      cosF.allocate(dimensions);
      cosFK.allocate(dimensions);
      HostDArray<cudaReal> cosF_h;
      cosF_h.allocate(meshSize);
      PolymerSpecies<cudaReal> const & polymer = system.mixture().polymerSpecies(0);
      for (int k = 0; k < meshSize; k++){
         cosF_h[k] = A * std::cos(2 * M_PI * k * f / meshSize);
         
         // Apply to each monomer
         if (PolymerModel::isBead()) {
            cosF_h[k] /= polymer.nBead();
         } else {
            cosF_h[k] /= polymer.length();
         }
      }
      cosF = cosF_h;
      
      // Convert to Fourier Space
      system.domain().fft().forwardTransform(cosF, cosFK);
      
      DArray< RField<1,CUT> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }      
      
      // Loop over monomer types add pressure perturbation
      DArray< RField<1,CUT> > const & w = system.w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::addVV(w2[i], w[i], cosF);
      }
      system.w().setRGrid(w2);
      
      system.compute();
      
      // Incompressibility error
      RField<1,CUT> error;
      error.allocate(dimensions);
       
      // Initialize resid to c field of species 0 minus 1
      VecOp::subVS(error, system.c().rgrid(0), 1.0);
      
      // Add other c fields to get SCF residual vector elements
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(error, system.c().rgrid(i));
      }
      
      // Intra analytical
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CUT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CUT> intra_(system);
      intra_.computeOmegaTotal(intraCorrelationK);
      
      // Compute analytical dphi using Intra
      RField<1,CUT> analyticalError;
      RFieldDft<1,CUT> analyticalErrorK;
      analyticalError.allocate(dimensions);
      analyticalErrorK.allocate(dimensions);
      VecOp::mulVV(analyticalErrorK, cosFK, intraCorrelationK);
      VecOp::mulEqS(analyticalErrorK, -1 * vMonomer);
      
      system.domain().fft().inverseTransformUnsafe(analyticalErrorK, analyticalError);
      
      RFieldComparison<1,CUT> comparison;
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
      Rp::System<1,CUT> system; 
      initSystem(system, "in/param_system_1D_diblcok_thread");
      system.w().readRGrid("in/w_diblock_homogenous.rf");
      
      IntVec<1> const & dimensions = system.domain().mesh().dimensions();
     
      // The intracorrelation function of conformational diblock 
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CUT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CUT> intra(system);
      intra.computeOmegaTotal(intraCorrelationK);
      
      // The intracorrelation function of conformational homo
      Rp::System<1,CUT> systemHomo; 
      initSystem(systemHomo, "in/param_system_1D_homo_thread");
      systemHomo.w().readRGrid("in/w_homo_homogenous.rf");
      RField<1,CUT> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CUT> intraHomo(systemHomo);
      intraHomo.computeOmegaTotal(intraCorrelationKHomo);
      
      RFieldComparison<1,CUT> comparison;
      comparison.compare(intraCorrelationK, intraCorrelationKHomo);
      TEST_ASSERT(comparison.maxDiff() < 1e-5);
      
   }
   
    void testIntraHomoBead()
   {
      // Compare the intracorrelation function of homopolymer and conformational diblock 
      printMethod(TEST_FUNC);
      
      openLogFile("out/testIntraHomoBead.log");
      Rp::System<1,CUT> system; 
      initSystem(system, "in/param_system_1D_diblcok_bead");
      system.w().readRGrid("in/w_diblock_homogenous.rf");
      
      IntVec<1> const & dimensions = system.domain().mesh().dimensions();
     
      // The intracorrelation function of conformational diblock 
      IntVec<1> kMeshDimensions;
      kMeshDimensions[0] = dimensions[0]/2 + 1;
      RField<1,CUT> intraCorrelationK;
      intraCorrelationK.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CUT> intra(system);
      intra.computeOmegaTotal(intraCorrelationK);
      
      // The intracorrelation function of conformational homo
      Rp::System<1,CUT> systemHomo; 
      initSystem(systemHomo, "in/param_system_1D_homo_bead");
      systemHomo.w().readRGrid("in/w_homo_homogenous.rf");
      RField<1,CUT> intraCorrelationKHomo;
      intraCorrelationKHomo.allocate(kMeshDimensions);
      Rp::IntraCorrelation<1,CUT> intraHomo(systemHomo);
      intraHomo.computeOmegaTotal(intraCorrelationKHomo);
      
      RFieldComparison<1,CUT> comparison;
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
