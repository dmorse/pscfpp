#ifndef RPC_COMPRESSOR_TEST_H
#define RPC_COMPRESSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rpc/fts/compressor/AmCompressor.h>
#include <rpc/fts/compressor/LrCompressor.h>
#include <rpc/fts/compressor/LrAmCompressor.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <prdc/field/cpu/RFieldComparison.h>

#include <pscf/mesh/Mesh.h>

#include <util/misc/FileMaster.h> 
#include <util/random/Random.h> 
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CompressorTest : public LogFileUnitTest
{

   Rp::System<3, CppTp<3> > system;
   
public:
   
   void setUp()
   {  setVerbose(0); }
   
   void tearDown()
   {  setVerbose(0); }
   
   template <int D>
   void initSystem(Rp::System<D, CppTp<D> >& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

   }
   
   template <int D>
   void randomStep(Rp::System<D, CppTp<D> >& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3, CppTp<3> > > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      
      DArray< RField<3, CppTp<3> > > const & w = system.w().rgrid();
      double stepSize = 1e-1;
      Random random;
      random.setSeed(0);
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            double r = random.uniform(-stepSize,stepSize);
            w2[i][k] =  w[i][k]+ r;
         }
      }
      system.w().setRGrid(w2);
   }
   
   template <int D>
   void addPressureField(Rp::System<D, CppTp<D> >& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3, CppTp<3> > > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      
      DArray< RField<3, CppTp<3> > > const & w = system.w().rgrid();
      double stepSize = 1e-1;
      Random random;
      random.setSeed(0);
      for (int k = 0; k < meshSize; k++){
         double r = random.uniform(-stepSize,stepSize);
         for (int i = 0; i < nMonomer; i++){
            w2[i][k] =  w[i][k]+ r;
         }
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
                       Rp::System<3, CppTp<3> >& system, 
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
      DArray< RField<3, CppTp<3> > > w0;
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
      TEST_ASSERT(sqrt(totalError)/sqrt(meshSize) < 1.0E-8);
      
      // Reset back to input chemical potential fields
      system.w().setRGrid(w0);
      
      // Apply pressure field
      addPressureField(system);
      compressor.compress();
      DArray< RField<3, CppTp<3> > > w1;
      w1.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w1[i].allocate(dimensions);
      }
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w1[i][j] = system.w().rgrid(i)[j];
         }
      }
      RFieldComparison<3, CppTp<3> > comparison;
      comparison.compare(w0, w1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-2);
   }
   
   
   void testAmCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3, CppTp<3> > system;
      Rp::AmCompressor<3, CppTp<3> > amCompressor(system);
      testCompressor(amCompressor, system, 
                     "in/param_AmCompressor",
                     "out/testAmCompressor.log");
   }
   
   void testLrCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3, CppTp<3> > system;
      Rp::LrCompressor<3, CppTp<3> > lrCompressor(system);
      testCompressor(lrCompressor,  system, 
                     "in/param_LrCompressor", 
                     "out/testLrCompressor.log");
   }

   #if 0
   void testLrAmPreCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3, CppTp<3> > system;
      LrAmPreCompressor<3> lrAmPreCompressor(system);
      testCompressor(lrAmPreCompressor, system, 
                     "in/param_LrAmPreCompressor",
                     "out/testLrAmPreCompressor.log");
   }
   #endif
   
   void testLrAmCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3, CppTp<3> > system;
      Rp::LrAmCompressor<3, CppTp<3> > lrAmCompressor(system);
      testCompressor(lrAmCompressor, system, 
                     "in/param_LrAmCompressor",
                     "out/testLrAmCompressor.log");
   }
   
};

TEST_BEGIN(CompressorTest)
TEST_ADD(CompressorTest, testAmCompressor)
TEST_ADD(CompressorTest, testLrCompressor)
//TEST_ADD(CompressorTest, testLrAmPreCompressor)
TEST_ADD(CompressorTest, testLrAmCompressor)
TEST_END(CompressorTest)

#endif
