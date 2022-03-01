#ifndef PSPG_TEST_CUDA_TEST_COMPOSITE_H
#define PSPG_TEST_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaMemTest.h"
#include "CudaResourceTest.h"

TEST_COMPOSITE_BEGIN(CudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaMemTest);
TEST_COMPOSITE_ADD_UNIT(CudaResourceTest);
TEST_COMPOSITE_END

#endif
