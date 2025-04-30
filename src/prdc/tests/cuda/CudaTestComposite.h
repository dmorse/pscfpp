#ifndef PRDC_CUDA_TEST_COMPOSITE_H
#define PRDC_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaReduceTest.h"
#include "CudaVecOpTest.h"
#include "CudaFieldTest.h"
#include "CudaFieldComparisonTest.h"
#include "CudaFieldTest.h"
#include "CudaFftTest.h"
#include "CudaComplexTest.h"
#include "WaveListTest.h"

TEST_COMPOSITE_BEGIN(CudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaReduceTest);
TEST_COMPOSITE_ADD_UNIT(CudaVecOpTest);
TEST_COMPOSITE_ADD_UNIT(CudaFieldTest);
TEST_COMPOSITE_ADD_UNIT(CudaFieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(CudaFftTest);
TEST_COMPOSITE_ADD_UNIT(CudaComplexTest);
TEST_COMPOSITE_ADD_UNIT(WaveListTest);
TEST_COMPOSITE_END

#endif
