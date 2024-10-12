#ifndef PRDC_CUDA_FIELD_TEST_H
#define PRDC_CUDA_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/Field.h>
#include <prdc/cuda/Field.tpp>
#include <prdc/cpu/Field.h>

using namespace Util;
using namespace Pscf; 
using namespace Pscf::Prdc; 

class CudaFieldTest : public UnitTest 
{

private:

   const static int capacity = 3;

public:

   void setUp(){}

   void tearDown(){}

   void testConstructor();
   void testAllocate();
   void testAssignRoundTrip1();

};


void CudaFieldTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      Cuda::Field<double> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void CudaFieldTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      Cuda::Field<double> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
   }
} 

void CudaFieldTest::testAssignRoundTrip1()
{
   printMethod(TEST_FUNC);
   {

      // Initialize host (cpu) data in field vh1
      Cpu::Field<double> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i] = (i+1)*10.0 ;
      }
  
      // Copy host field vh1 to device field vd
      Cuda::Field<double> vd;
      //vd.allocate(capacity);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());
       
      // Copy device field vd to host host field vh2
      Cpu::Field<double> vh2;
      vh2 = vd;

      TEST_ASSERT(vh2.capacity() == vh1.capacity());
      TEST_ASSERT(vh2[0] == 10.0);
      TEST_ASSERT(vh2[2] == 30.0);
      for (int i=0; i < capacity; i++ ) {
         TEST_ASSERT(eq(vh2[i], (i+1)*10.0)) ;
      }
   }
} 

TEST_BEGIN(CudaFieldTest)
TEST_ADD(CudaFieldTest, testConstructor)
TEST_ADD(CudaFieldTest, testAllocate)
TEST_ADD(CudaFieldTest, testAssignRoundTrip1)
TEST_END(CudaFieldTest)
#endif
