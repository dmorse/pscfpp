#ifndef RPG_COMPRESSOR_TEST_H
#define RPG_COMPRESSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/fts/compressor/AmCompressor.h>
#include <rpg/fts/compressor/LrCompressor.h>
#include <rpg/fts/compressor/LrAmPreCompressor.h>
#include <rpg/fts/compressor/LrAmCompressor.h>

#include <prdc/cuda/RFieldComparison.h>

#include <pscf/cuda/CudaRandom.h> 
#include <pscf/cuda/GpuResources.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/random/Random.h> 

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class CompressorTest : public LogFileUnitTest
{

   System<3> system;
   
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
   
   template <int D>
   void randomStep(System<D>& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      RField<D> randomField;
      randomField.allocate(dimensions);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      double stepSize = 1e-1;
      CudaRandom cudaRandom;
      cudaRandom.setSeed(0);
      DArray< RField<3> > const & w = system.w().rgrid();
      
      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++){

         // Generate random numbers between 0.0 and 1.0 from uniform distribution
         cudaRandom.uniform(randomField.cArray(), meshSize);

         // Generate random numbers between [-stepSize_,stepSize_]
         mcftsScale<<<nBlocks, nThreads>>>(randomField.cArray(), stepSize, meshSize);

         // Change the w field configuration
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>(w[i].cArray(), randomField.cArray(),
                                                   w2[i].cArray(), meshSize);

      }
      
      // set system r grid
      system.setWRGrid(w2);

   }
   
   template <int D>
   void addPressureField(System<D>& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      DArray< RField<3> > const & w = system.w().rgrid();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      RField<D> randomField;
      randomField.allocate(dimensions);
      
      CudaRandom cudaRandom;
      cudaRandom.setSeed(0);
      cudaRandom.uniform(randomField.cArray(), meshSize);
      double stepSize = 1e-1;
      mcftsScale<<<nBlocks, nThreads>>>(randomField.cArray(), stepSize, meshSize);
      
      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++){
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>(w[i].cArray(), randomField.cArray(),
                                                   w2[i].cArray(), meshSize);
      }
      system.setWRGrid(w2);
      
   }
   
   template <typename Compressor>
   void initCompressor(Compressor& compressor, std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      compressor.readParam(in);
      in.close();
   }
  
   /*
   * Generic test function template.
   */ 
   template <typename Compressor>
   void testCompressor(Compressor& compressor, 
                       System<3>& system, 
                       std::string infilename, 
                       char const * outfilename)
   {
      openLogFile(outfilename);
      
      initSystem(system, "in/param_system_disordered");
      initCompressor(compressor, infilename);
      system.readWRGrid("in/w_dis.rf");
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<3> const & dimensions = system.domain().mesh().dimensions();
      
      // Store value of input chemical potential fields
      DArray< RField<3> > w0;
      w0.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0[i].allocate(dimensions);
      }
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      DArray<RField<3>> const * currSys = &system.w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w0[i].cArray(), 
                                          (*currSys)[i].cArray(), meshSize);
         
      }
      
      // Apply a random step to check the incompressibility constraint
      randomStep(system);
      
      compressor.compress();
      
      // Compute incompressible error
      RField<3> error;
      error.allocate(dimensions);
      assignUniformReal<<<nBlocks, nThreads>>>(error.cArray(), -1.0, meshSize);
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (error.cArray(), system.c().rgrid(i).cArray(), meshSize);
      }
      double product = (double)gpuInnerProduct(error.cArray(), error.cArray(), meshSize);
      
      TEST_ASSERT(sqrt(product)/sqrt(meshSize) < 1.0E-8);
      
      // Reset back to input chemical potential fields
      system.setWRGrid(w0);
      
      // Apply pressure field
      addPressureField(system);
      
      compressor.compress();
      DArray< RField<3> > w1;
      w1.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w1[i].allocate(dimensions);
      }
      
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w1[i].cArray(), 
                                          (*currSys)[i].cArray(), meshSize);
      }
      
      RFieldComparison<3> comparison;
      comparison.compare(w0, w1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-2);
      
   }
   
   
   void testAmCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      AmCompressor<3> amCompressor(system);
      testCompressor(amCompressor, system, 
                     "in/param_AmCompressor",
                     "out/testAmCompressor.log");
   }
   
   void testLrCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      LrCompressor<3> lrCompressor(system);
      testCompressor(lrCompressor,  system, 
                     "in/param_LrCompressor", 
                     "out/testLrCompressor.log");
   }

   void testLrAmPreCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      LrAmPreCompressor<3> lrAmPreCompressor(system);
      testCompressor(lrAmPreCompressor, system, 
                     "in/param_LrAmPreCompressor",
                     "out/testLrAmPreCompressor.log");
   }
   
   void testLrAmCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      LrAmCompressor<3> lrAmCompressor(system);
      testCompressor(lrAmCompressor, system, 
                     "in/param_LrAmCompressor",
                     "out/testLrAmCompressor.log");
   }
   
};

TEST_BEGIN(CompressorTest)
TEST_ADD(CompressorTest, testAmCompressor)
TEST_ADD(CompressorTest, testLrCompressor)
TEST_ADD(CompressorTest, testLrAmPreCompressor)
TEST_ADD(CompressorTest, testLrAmCompressor)
TEST_END(CompressorTest)

#endif
