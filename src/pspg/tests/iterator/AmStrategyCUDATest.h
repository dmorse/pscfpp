#ifndef PSPG_AM_STRATEGY_CUDA_TEST_H
#define PSPG_AM_STRATEGY_CUDA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/iterator/AmStrategyCUDA.h>

#include <pspg/field/RDField.h>
#include <util/math/Constants.h>
#include <pspg/math/GpuResources.h>

#include <fstream>
#include <iomanip>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class AmStrategyCUDATest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}
   

};

TEST_BEGIN(AmStrategyCUDATest)
TEST_END(AmStrategyCUDATest)
#endif
