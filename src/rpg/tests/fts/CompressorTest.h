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
#include <rpg/fts/compressor/LrAmCompressor.h>
#include <rpg/fts/VecOpFts.h>

#include <prdc/cuda/RFieldComparison.h>
#include <prdc/cuda/resources.h>

#include <pscf/cuda/CudaRandom.h> 

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
      
      double stepSize = 1e-1;
      CudaRandom cudaRandom;
      cudaRandom.setSeed(0);
      DArray< RField<3> > const & w = system.w().rgrid();
      
      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++){

         // Generate random numbers between 0.0 and 1.0 from uniform distribution
         cudaRandom.uniform(randomField);

         // Generate random numbers between [-stepSize_,stepSize_]
         VecOpFts::mcftsScale(randomField, stepSize);

         // Change the w field configuration
         VecOp::addVV(w2[i], w[i], randomField);

      }
      
      // set system r grid
      system.w().setRGrid(w2);

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
      
      RField<D> randomField;
      randomField.allocate(dimensions);
      
      CudaRandom cudaRandom;
      cudaRandom.setSeed(0);
      cudaRandom.uniform(randomField);
      double stepSize = 1e-1;
      VecOpFts::mcftsScale(randomField, stepSize);
      
      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(w2[i], w[i], randomField);
      }
      system.w().setRGrid(w2);
      
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
      system.w().readRGrid("in/w_dis.rf");
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<3> const & dimensions = system.domain().mesh().dimensions();
      
      // Store value of input chemical potential fields
      DArray< RField<3> > w0;
      w0.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0[i].allocate(dimensions);
      }
      
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0[i], system.w().rgrid(i));
         
      }
      
      // Apply a random step to check the incompressibility constraint
      randomStep(system);
      
      compressor.compress();
      
      // Compute incompressible error
      RField<3> error;
      error.allocate(dimensions);
      VecOp::eqS(error, -1.0);
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addEqV(error, system.c().rgrid(i));
      }
      double product = Reduce::innerProduct(error, error); 
      
      TEST_ASSERT(sqrt(product)/sqrt(meshSize) < 1.0E-8);
      
      // Reset back to input chemical potential fields
      system.w().setRGrid(w0);
      
      // Apply pressure field
      addPressureField(system);
      
      compressor.compress();
      DArray< RField<3> > w1;
      w1.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w1[i].allocate(dimensions);
      }
      
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w1[i], system.w().rgrid(i));
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
TEST_ADD(CompressorTest, testLrAmCompressor)
TEST_END(CompressorTest)

#endif
