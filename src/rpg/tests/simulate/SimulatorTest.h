#ifndef RPG_SIMULATOR_TEST_H
#define RPG_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/fts/Simulator.h>
#include <prdc/cuda/RField.h> 
#include <util/containers/DArray.h>  

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class SimulatorTest : public LogFileUnitTest
{

   System<3> system;

public:


   SimulatorTest()
    : system()
   {}

   void setUp()
   {  setVerbose(0); }

   void initSystem(std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testSystem.log");
      ParamComponent::setEcho(true);

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

   }

   void testAnalyzeChi()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param1_simulator");
      Simulator<3> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      double chi = system.interaction().chi(0,1);
      TEST_ASSERT( fabs(system.interaction().chi(0,0)) < 1.0E-8);
      TEST_ASSERT( fabs(system.interaction().chi(1,1)) < 1.0E-8);

      DArray<double> vals = simulator.chiEvals();
      TEST_ASSERT( fabs((vals[0] - simulator.chiEval(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((vals[1] - simulator.chiEval(1))/chi) < 1.0E-8);
      TEST_ASSERT(fabs((vals[0] + chi)/chi) < 1.0E-8);
      TEST_ASSERT(fabs(vals[1]/chi) < 1.0E-8);

      DMatrix<double> vecs = simulator.chiEvecs();
      TEST_ASSERT(fabs(vecs(0,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(0,1) + 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,1) - 1.0) < 1.0E-8);

      DArray<double>  sc = simulator.sc();
      TEST_ASSERT( fabs((sc[0] - simulator.sc(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((sc[1] - simulator.sc(1))/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(0)/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(1)/chi - 0.5) < 1.0E-8);

      #if 0
      std::cout << std::endl;
      std::cout << "vals  = " << vals[0] << "  " << vals[1] << std::endl;
      std::cout << "vec0  = " << vecs(0,0) << "  " << vecs(0,1) << std::endl;
      std::cout << "vec1  = " << vecs(1,0) << "  " << vecs(1,1) << std::endl;
      #endif

   }

   void testSaddlePointField()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param1_simulator");
      Simulator<3> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      system.readWRGrid("in/w_gyr.rf");
      DArray< RField<3> > const & w = system.w().rgrid();

      system.compute();
      DArray< RField<3> > const & c = system.c().rgrid();

      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();

      simulator.computeWc();
      DArray< RField<3> > const & wc = simulator.wc();

      simulator.computeCc();
      DArray< RField<3> > const & cc = simulator.cc();

      simulator.computeDc();
      DArray< RField<3> > const & dc = simulator.dc();

      // Check allocation and capacities
      TEST_ASSERT(c.capacity() == nMonomer);
      TEST_ASSERT(w.capacity() == nMonomer);
      TEST_ASSERT(wc.capacity() == nMonomer);
      TEST_ASSERT(cc.capacity() == nMonomer);
      int i;
      for (i=0; i < nMonomer; ++i) {
         TEST_ASSERT(c[i].capacity() == meshSize);
         TEST_ASSERT(w[i].capacity() == meshSize);
         TEST_ASSERT(wc[i].capacity() == meshSize);
         TEST_ASSERT(cc[i].capacity() == meshSize);
      }
      TEST_ASSERT(dc.capacity() == nMonomer - 1);
      for (i=0; i < nMonomer - 1; ++i) {
         TEST_ASSERT(dc[i].capacity() == meshSize);
      }

      // Test wc field
      RField<3> wcTest1;
      wcTest1.allocate(meshSize);
      RField<3> wcTest2;
      wcTest2.allocate(meshSize);
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      /// TEST_ASSERT(fabs( w[0][i] - wc[0][i] - wc[1][i] ) < 1.0E-6);
      pointWiseBinarySubtract<<<nBlocks, nThreads>>>
         (w[0].cField(), wc[0].cField(), wcTest1.cField(), meshSize);
      pointWiseSubtract<<<nBlocks, nThreads>>>
         (wcTest1.cField(), wc[1].cField(), meshSize);
      TEST_ASSERT((double)gpuMaxAbs(wcTest1.cField(), meshSize)< 1.0E-6);
      
      /// TEST_ASSERT(fabs( w[0][i] - w[1][i] - 2.0*wc[0][i] ) < 1.0E-6);
      pointWiseBinarySubtract<<<nBlocks, nThreads>>>
         (w[0].cField(), w[1].cField(), wcTest2.cField(), meshSize);
      pointWiseAddScale<<<nBlocks, nThreads>>>
         (wcTest2.cField(), wc[0].cField(), -2.0, meshSize);
      TEST_ASSERT((double)gpuMaxAbs(wcTest2.cField(), meshSize)< 1.0E-6);


      // Test cc field
      RField<3> ccTest;
      ccTest.allocate(meshSize);
      ///TEST_ASSERT(fabs( c[0][i] - c[1][i] - cc[0][i] ) < 1.0E-6);
      pointWiseBinarySubtract<<<nBlocks, nThreads>>>
         (c[0].cField(), c[1].cField(), ccTest.cField(), meshSize);
      pointWiseSubtract<<<nBlocks, nThreads>>>
         (ccTest.cField(), cc[0].cField(), meshSize);
      TEST_ASSERT((double)gpuMaxAbs(ccTest.cField(), meshSize)< 1.0E-6);
      
      // Test dc field
      // diff = fabs( dc[0][i] );
      TEST_ASSERT((double)gpuMaxAbs(dc[0].cField(), meshSize) < 1.0E-8);

      #if 0
      std::cout << std::endl;
      std::cout << "Maximum derivative " << diff << std::endl;
      std::cout << std::endl;
      std::cout << system.domain().mesh().dimensions() << std::endl;
      std::cout << system.domain().mesh().size() << std::endl;
      #endif

      double volume = system.domain().unitCell().volume();
      double vMonomer = system.mixture().vMonomer();
      double ratio = volume/vMonomer;

      // SCFT free energy for converged solution
      system.computeFreeEnergy();
      double fHelmholtz = system.fHelmholtz();

      // FTS Hamiltonian at saddle-point
      simulator.computeHamiltonian();
      double hamiltonian = simulator.hamiltonian();

      // Compare FTS Hamiltonian to SCFT free energy
      double diff;
      diff = fabs((hamiltonian - ratio*fHelmholtz)/hamiltonian);
      TEST_ASSERT( diff < 1.0E-8);

      #if 0
      std::cout << "Hamiltonian different (fractional) = " 
                << diff << std::endl;

      std::cout << "fHelmholtz = " << fHelmholtz << "  " 
                << ratio*fHelmholtz  << std::endl;

      std::cout << "Hamiltonian = " << hamiltonian/ratio << "  " 
                << hamiltonian << std::endl;
      #endif
     
   }
};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testAnalyzeChi)
TEST_ADD(SimulatorTest, testSaddlePointField)
TEST_END(SimulatorTest)

#endif
