#ifndef PSPC_TEST_CUDA_TEST_COMPOSITE_H
#define PSPC_TEST_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaMemTest.h"

TEST_COMPOSITE_BEGIN(CudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaMemTest);
TEST_COMPOSITE_END

#endif
