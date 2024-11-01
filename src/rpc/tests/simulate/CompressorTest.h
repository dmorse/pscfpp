#ifndef RPC_COMPRESSOR_TEST_H
#define RPC_COMPRESSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>
#include <rpc/fts/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/compressor/AmCompressor.h>
#include <rpc/fts/compressor/LrCompressor.h>
#include <rpc/fts/compressor/LrAmCompressor.h>

#include <prdc/cpu/RFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/random/Random.h> 

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

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
      
      DArray< RField<3> > const & w = system.w().rgrid();
      double stepSize = 1e-2;
      Random random;
      random.setSeed(12345);
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            double r = random.uniform(-stepSize,stepSize);
            w2[i][k] =  w[i][k]+ r;
         }
      }
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
      double stepSize = 1e-2;
      Random random;
      random.setSeed(12345);
      for (int k = 0; k < meshSize; k++){
         double r = random.uniform(-stepSize,stepSize);
         for (int i = 0; i < nMonomer; i++){
            w2[i][k] =  w[i][k]+ r;
         }
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
   
   template <typename Compressor>
   void testCompressor(Compressor& compressor, System<3>& system, std::string infilename, char const * outfilename)
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
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < meshSize; ++j){
            w0[i][j] = system.w().rgrid(i)[j];
         }
      }
      
      // Apply a random step to check the incompressibility constraint
      randomStep(system);
      compressor.compress();
      double totalError = 0.0;
      for (int i = 0; i<  meshSize; i++){
         double error = -1.0;
         for (int j = 0; j <nMonomer ; j++){
            error +=  system.c().rgrid(j)[i];
         }
         totalError += error*error;
      }
      TEST_ASSERT(sqrt(totalError)/sqrt(meshSize) < 1.0E-4);
      
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
         for (int j = 0; j< meshSize; ++j){
            w1[i][j] = system.w().rgrid(i)[j];
         }
      }
      RFieldComparison<3> comparison;
      comparison.compare(w0, w1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-4);
   }
   
   
   void testAmCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      AmCompressor<3> amCompressor(system);
      testCompressor(amCompressor, system, "in/param_AmCompressor","out/testAmCompressor.log");
   }
   
   void testLrAmCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      LrAmCompressor<3> lrAmCompressor(system);
      testCompressor(lrAmCompressor, system, "in/param_LrAmCompressor","out/testLrAmCompressor.log");
   }
   
   void testLrCompressor()
   {
      printMethod(TEST_FUNC);
      System<3> system;
      LrCompressor<3> lrCompressor(system);
      testCompressor(lrCompressor,  system, "in/param_LrCompressor", "out/testLrCompressor.log");
   }

};

TEST_BEGIN(CompressorTest)
TEST_ADD(CompressorTest, testAmCompressor)
TEST_ADD(CompressorTest, testLrAmCompressor)
TEST_ADD(CompressorTest, testLrCompressor)
TEST_END(CompressorTest)

#endif
