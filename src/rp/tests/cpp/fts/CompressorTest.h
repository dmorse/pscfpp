#ifndef RPC_COMPRESSOR_TEST_H
#define RPC_COMPRESSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/fts/compressor/AmCompressor.h>
#include <rp/fts/compressor/LrCompressor.h>
#include <rp/fts/compressor/LrAmCompressor.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/field/cpp/RFieldComparison.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>

#include <util/misc/FileMaster.h> 
#include <util/random/Random.h> 
#include <util/math/Constants.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CompressorTest : public LogFileUnitTest
{

   Rp::System<3,CPT> system;
   
public:
   
   void setUp()
   {  setVerbose(0); }
   
   void tearDown()
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
   
   template <int D>
   void randomStep(Rp::System<D,CPT>& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3,CPT> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      
      DArray< RField<3,CPT> > const & w = system.w().rgrid();
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
   void addPressureField(Rp::System<D,CPT>& system)
   {
      // Random change in pressure field
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      IntVec<D> const & dimensions = system.domain().mesh().dimensions();
      DArray< RField<3,CPT> > w2;
      w2.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w2[i].allocate(dimensions);
      }
      
      DArray< RField<3,CPT> > const & w = system.w().rgrid();
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
                       Rp::System<3,CPT>& system, 
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
      DArray< RField<3,CPT> > w0;
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
      DArray< RField<3,CPT> > w1;
      w1.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w1[i].allocate(dimensions);
      }
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w1[i][j] = system.w().rgrid(i)[j];
         }
      }
      RFieldComparison<3,CPT> comparison;
      comparison.compare(w0, w1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-2);
   }
   
   
   void testAmCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      Rp::AmCompressor<3,CPT> amCompressor(system);
      testCompressor(amCompressor, system, 
                     "in/param_AmCompressor",
                     "out/testAmCompressor.log");
   }
   
   void testLrCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      Rp::LrCompressor<3,CPT> lrCompressor(system);
      testCompressor(lrCompressor,  system, 
                     "in/param_LrCompressor", 
                     "out/testLrCompressor.log");
   }

   #if 0
   void testLrAmPreCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      LrAmPreCompressor<3> lrAmPreCompressor(system);
      testCompressor(lrAmPreCompressor, system, 
                     "in/param_LrAmPreCompressor",
                     "out/testLrAmPreCompressor.log");
   }
   #endif
   
   void testLrAmCompressor()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      Rp::LrAmCompressor<3,CPT> lrAmCompressor(system);
      testCompressor(lrAmCompressor, system, 
                     "in/param_LrAmCompressor",
                     "out/testLrAmCompressor.log");
   }

   /*
   * Commands after a COMPRESS command that failed must not be executed.
   */
   void testCompressCommandFailure()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testCompressCommandFailure.log");

      // Diblock melt, with an AmCompressor that cannot converge
      Rp::System<1,CPT> system1;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.02 }"
         " Interaction{ chi( 1 0 12.0 ) }"
         " Domain{ mesh 16 lattice lamellar }"
         " BdSimulator{ seed 1234 LMBdStep{ mobility 0.01 }"
         " AmCompressor{ epsilon 1.0e-12 maxItr 1 maxHist 10 } } }\n");
      system1.readParam(param);
      TEST_ASSERT(system1.hasSimulator());
      FSArray<double, 6> parameters;
      parameters.append(1.0);
      system1.setUnitCell(parameters);
      DArray< RField<1,CPT> > w;
      w.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(system1.domain().mesh().dimensions());
         for (int k = 0; k < 16; ++k) {
            w[i][k] = (i == 0) ? cos(2.0*Constants::Pi*k/16.0) : 0.0;
         }
      }
      system1.w().setRGrid(w);

      std::istringstream in("COMPRESS\nSET_UNIT_CELL lamellar 1.5\n"
                            "FINISH\n");
      TEST_ASSERT(system1.readCommands(in) == 1);
      TEST_ASSERT(std::abs(system1.domain().unitCell().parameter(0) - 1.0)
                  < 1.0E-12);
   }

   /*
   * Set the w fields of a 1D system to w_0 = a + b*cos(2 pi n x/L) and
   * w_i = a for i > 0.
   */
   void setW1D(Rp::System<1,CPT>& system, double a, double b, int n)
   {
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      DArray< RField<1,CPT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < meshSize; ++k) {
            w[i][k] = a;
            if (i == 0) w[i][k] += b*std::cos(2.0*M_PI*n*k/meshSize);
         }
      }
      system.w().setRGrid(w);
   }

   /*
   * Root-mean-square deviation of the total concentration from 1.
   */
   double rmsResidual(Rp::System<1,CPT>& system)
   {
      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();
      double sumSq = 0.0;
      for (int k = 0; k < meshSize; ++k) {
         double error = -1.0;
         for (int i = 0; i < nMonomer; ++i) {
            error += system.c().rgrid(i)[k];
         }
         sumSq += error*error;
      }
      return std::sqrt(sumSq/double(meshSize));
   }

   /*
   * COMPRESS with an LrCompressor after the unit cell has changed.
   */
   void testLrCompressorCellChange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLrCompressorCellChange.log");

      // Diblock melt, with an LrCompressor
      Rp::System<1,CPT> system;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.01 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar }"
         " BdSimulator{ seed 1234 LMBdStep{ mobility 0.01 }"
         " LrCompressor{ epsilon 1.0e-9 maxItr 100 } } }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(1.0);
      system.setUnitCell(parameters);
      setW1D(system, 0.0, 1.0, 1);

      // Compress, triple the cell length, and compress again
      std::istringstream in("COMPRESS\nSET_UNIT_CELL lamellar 3.0\n"
                            "COMPRESS\nFINISH\n");
      system.readCommands(in);
      TEST_ASSERT(system.c().hasData());
      TEST_ASSERT(rmsResidual(system) < 1.0E-8);
   }

   /*
   * Blend of two open homopolymers A and B, in fields w_A = a + 0.2*cos
   * and w_B = a, so that the total volume fraction is about exp(-a).
   */
   void initOpenBlend(Rp::System<1,CPT>& system, double a)
   {
      // The volume fraction of each species is 0.5*exp(-w) in a
      // homogeneous field w
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 2"
         " Polymer{ type linear nBlock 1 blocks[ 0 1.0 ] mu -0.6931 }"
         " Polymer{ type linear nBlock 1 blocks[ 1 1.0 ] mu -0.6931 }"
         " vMonomer 0.01 ds 0.02 }"
         " Interaction{ chi( 1 0 2.0 ) }"
         " Domain{ mesh 16 lattice lamellar } }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      system.setUnitCell(parameters);
      setW1D(system, a, 0.2, 1);
   }

   /*
   * LR compressors in a blend of open species.
   */
   void testLrCompressorOpen()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLrCompressorOpen.log");

      // LrCompressor, from a total volume fraction exp(-1): the
      // correlation must follow the volume fractions at each step
      Rp::System<1,CPT> system1;
      initOpenBlend(system1, 1.0);
      Rp::LrCompressor<1,CPT> lrCompressor(system1);
      std::istringstream in1("LrCompressor{ epsilon 1.0e-9 maxItr 100 }\n");
      lrCompressor.readParam(in1);
      TEST_ASSERT(lrCompressor.compress() == 0);
      TEST_ASSERT(rmsResidual(system1) < 1.0E-8);

      // LrAmCompressor, from total volume fractions exp(-0.3) and
      // exp(-1)
      double a[2] = {0.3, 1.0};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system2;
         initOpenBlend(system2, a[i]);
         Rp::LrAmCompressor<1,CPT> lrAmCompressor(system2);
         std::istringstream in2("LrAmCompressor{ epsilon 1.0e-9"
                                " maxItr 100 maxHist 10 }\n");
         lrAmCompressor.readParam(in2);
         TEST_ASSERT(lrAmCompressor.compress() == 0);
         TEST_ASSERT(rmsResidual(system2) < 1.0E-8);
         if (i == 1) {
            // Second compression by the same compressor, from exp(1)
            setW1D(system2, -1.0, 0.2, 1);
            TEST_ASSERT(lrAmCompressor.compress() == 0);
            TEST_ASSERT(rmsResidual(system2) < 1.0E-8);
         }
      }
   }
   
};

TEST_BEGIN(CompressorTest)
TEST_ADD(CompressorTest, testAmCompressor)
TEST_ADD(CompressorTest, testLrCompressor)
//TEST_ADD(CompressorTest, testLrAmPreCompressor)
TEST_ADD(CompressorTest, testLrAmCompressor)
TEST_ADD(CompressorTest, testCompressCommandFailure)
TEST_ADD(CompressorTest, testLrCompressorCellChange)
TEST_ADD(CompressorTest, testLrCompressorOpen)
TEST_END(CompressorTest)

#endif
