#ifndef PRDC_CPU_FIELD_BASIS_CONVERTER__TEST_H
#define PRDC_CPU_FIELD_BASIS_CONVERTER__TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/FieldBasisConverter.h>
#include <prdc/cpu/RField.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

using namespace Util;
using namespace Pscf::Prdc::Cpu;

class CpuFieldBasisConverterTest : public UnitTest 
{

private:

   const static int nMonomer = 3;
   const static int meshSize = 3;

   IntVec<1> meshDimensions;
   DMatrix<double> basis;

   typedef double Data;

   DArray< RField<1> > ref;
   DArray< RField<1> > in;
   DArray< RField<1> > out;

public:

   void setUp() 
   {
      meshDimensions[0] = meshSize;
      basis.allocate(nMonomer, nMonomer);
     
      double prefactor = sqrt(3.0/6.0);
      basis(0,0) =  1.0*prefactor;
      basis(0,1) = -2.0*prefactor;
      basis(0,2) =  1.0*prefactor;

      prefactor = sqrt(3.0/2.0);
      basis(1,0) =  1.0*prefactor;
      basis(1,1) =  0.0*prefactor;
      basis(1,2) = -1.0*prefactor;

      basis(2,0) =  1.0;
      basis(2,1) =  1.0;
      basis(2,2) =  1.0;
      
   }

   void tearDown() {}

   void allocate() 
   {
      in.allocate(nMonomer);
      out.allocate(nMonomer);
      ref.allocate(nMonomer);
      int i, k;
      for (i = 0; i < nMonomer; ++i) {
         in[i].allocate(meshDimensions);
         out[i].allocate(meshDimensions);
         ref[i].allocate(meshDimensions);
      }
      in[0][0] =  0.7;
      in[0][1] = -0.3;
      in[0][2] = -0.2;
      in[1][0] = -0.2;
      in[1][1] = -0.9;
      in[1][2] =  0.1;
      in[2][0] =  0.8;
      in[2][1] = -0.3;
      in[2][2] =  0.5;
      for (i = 0; i < nMonomer; ++i) {
         for (k = 0; k < meshSize; ++k) {
            ref[i][k] = in[i][k];
            out[i][k] = 0.0;
         }
      }
   }

   void testDefaultConstructor();
   void testConstructor();
   void testRoundTrip();

};


void CpuFieldBasisConverterTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      FieldBasisConverter<1> converter;
   }
}

void CpuFieldBasisConverterTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      FieldBasisConverter<1> converter(basis);

      // Test Basis
      double normSq = 3.0;
      double error = converter.maxBasisError(normSq);
      //std::cout << "Max error = " << error << std::endl;
      TEST_ASSERT(error >= 0.0);
      TEST_ASSERT(fabs(error) < 1.0E-8);
   }
}

void CpuFieldBasisConverterTest::testRoundTrip()
{
   printMethod(TEST_FUNC);
   {
      FieldBasisConverter<1> converter(basis);
      allocate();

      // Round trip conversion 
      double prefactor = 1.0/double(nMonomer);
      converter.convertToBasis(in, out, prefactor);
      converter.convertFromBasis(out, in, 1.0);

      // Compare to reference array
      int i, j;
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < meshSize; ++j) {
            TEST_ASSERT(fabs(ref[i][j] - in[i][j]) < 1.0E-8);
         }
      }
   }

}

TEST_BEGIN(CpuFieldBasisConverterTest)
TEST_ADD(CpuFieldBasisConverterTest, testDefaultConstructor)
TEST_ADD(CpuFieldBasisConverterTest, testConstructor)
TEST_ADD(CpuFieldBasisConverterTest, testRoundTrip)
TEST_END(CpuFieldBasisConverterTest)

#endif
