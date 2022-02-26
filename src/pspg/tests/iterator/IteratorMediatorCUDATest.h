#ifndef PSPG_ITERATOR_MEDIATOR_CUDA_TEST_H
#define PSPG_ITERATOR_MEDIATOR_CUDA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/iterator/IteratorMediatorCUDA.h>

#include <pspg/field/RDField.h>
#include <util/math/Constants.h>
#include <pspg/math/GpuResources.h>

#include <fstream>
#include <iomanip>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class IteratorMediatorCUDATest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testNElements() 
   {
      printMethod(TEST_FUNC);

      System<3> system;
      setupSystem3D(system);

      IteratorMediatorCUDA<3> itermed(system);

      int n = itermed.nElements();

      int nCheck = system.mesh().size()*system.mixture().nMonomer();

      TEST_ASSERT(n==nCheck);
   }

   void testGetCurrent()
   {
      printMethod(TEST_FUNC);

      // Sample system
      System<3> system;
      setupSystem3D(system);
      IteratorMediatorCUDA<3> itermed(system);

      // Relevant parameters
      const int nMonomer = system.mixture().nMonomer();
      const int nMesh = system.mesh().size();
      const int n = itermed.nElements();

      // A field to store IteratorMediator output in
      FieldCUDA d_curr;
      d_curr.allocate(n);

      // Host arrays to store output in
      cudaReal* currCheck = new cudaReal[n];
      cudaReal* curr = new cudaReal[n];

      // Run IteratorMediator function. Note that this, per the original pspg convention,
      // modifies the field by subtracting out the average.
      itermed.getCurrent(d_curr);
      cudaMemcpy(curr, d_curr.cDField(), n*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      // Manually get current field
      for (int i = 0; i < nMonomer; i++) {
         cudaMemcpy(currCheck + i*nMesh, system.wFieldsRGrid()[i].cDField(), 
                              nMesh*sizeof(cudaReal), cudaMemcpyDeviceToHost);
      }

      // compare results
      for (int i = 0; i < n; i++) {
         // std::cout << "IterMed = " << curr[i] 
         //          << "    Test = " << currCheck[i] << std::endl;
         TEST_ASSERT(curr[i]==currCheck[i]);
      }
      
   }

   void testGetResidual()
   {
      printMethod(TEST_FUNC);

      // Sample system
      System<3> system;
      setupSystem3D(system);
      IteratorMediatorCUDA<3> itermed(system);

      // Relevant parameters
      const int nMonomer = system.mixture().nMonomer();
      const int nMesh = system.mesh().size();
      const int n = itermed.nElements();

      // A field to store IteratorMediator output in
      FieldCUDA d_resid;
      d_resid.allocate(n);

      // Host arrays to store output in
      cudaReal* resid = new cudaReal[n];

      // Run IteratorMediator function
      itermed.evaluate();
      itermed.getResidual(d_resid);
      cudaMemcpy(resid, d_resid.cDField(), n*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      // compare results
      bool passed = false;
      for (int i = 0; i < n; i++) {
         //std::cout << "resid = " << resid[i] << std::endl;
         if (fabs(resid[i]) > 1E-15) passed = true;
      }
      if (!passed)
         TEST_THROW("All residual elements were zero. This is not correct, so something went wrong.");

   }

   #if 0
   void testFindAverage()
   {
      printMethod(TEST_FUNC);

      // Sample system
      System<3> system;
      setupSystem3D(system);
      IteratorMediatorCUDA<3> itermed(system);

      // create host and device data arrays
      const int n = 128*128;
      cudaReal* data = new cudaReal[n];
      FieldCUDA d_data;
      d_data.allocate(n);

      // Create random data, store on host and device
      for (int i = 0; i < n; i++ ) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cDField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Compute average on host
      cudaReal avgTest = 0;
      for (int i = 0; i < n; i++) {
         avgTest += data[i];
      }
      avgTest/=n;

      // Compute average on device
      cudaReal avg;
      avg = itermed.findAverage(d_data.cDField(),n);
      
      std::cout << "Host Average = " << avgTest << std::endl;
      std::cout << "Device Average = " << avg << std::endl; 
      std::cout << "Difference = " << avg-avgTest << std::endl;

      TEST_ASSERT(avg==avgTest);
   }
   #endif

   void setupSystem3D(System<3> & system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());


      std::ifstream in;
      openInputFile("in/diblock/bcc/param.rigid", in);
      system.readParam(in);
      in.close();

      system.readWBasis("in/diblock/bcc/omega.in");
   }
   

};

TEST_BEGIN(IteratorMediatorCUDATest)
TEST_ADD(IteratorMediatorCUDATest, testNElements)
TEST_ADD(IteratorMediatorCUDATest, testGetCurrent)
TEST_ADD(IteratorMediatorCUDATest, testGetResidual)
// This test is currently disabled because the findAverage member should be private.
// TEST_ADD(IteratorMediatorCUDATest, testFindAverage)
TEST_END(IteratorMediatorCUDATest)

#endif
