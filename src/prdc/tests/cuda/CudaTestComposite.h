#ifndef PRDC_CUDA_TEST_COMPOSITE_H
#define PRDC_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaFieldTest.h"
#include "CudaFieldComparisonTest.h"
#include "CudaResourceTest.h"
#include "CudaFieldTest.h"
#include "CudaFftTest.h"

TEST_COMPOSITE_BEGIN(CudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaFieldTest);
TEST_COMPOSITE_ADD_UNIT(CudaFieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(CudaResourceTest);
TEST_COMPOSITE_ADD_UNIT(CudaFftTest);
TEST_COMPOSITE_END

#endif
