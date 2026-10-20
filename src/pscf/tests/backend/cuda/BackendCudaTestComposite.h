#ifndef PSCF_BACKEND_CUDA_TEST_COMPOSITE_H
#define PSCF_BACKEND_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaArrayTest.h"
#include "CudaVecRandomTest.h"
#include "CudaThreadGridTest.h"
#include "CudaComplexTest.h"
#include "CudaVecOpTest.h"
#include "CudaReduceTest.h"

TEST_COMPOSITE_BEGIN(BackendCudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaArrayTest);
TEST_COMPOSITE_ADD_UNIT(CudaVecRandomTest);
TEST_COMPOSITE_ADD_UNIT(CudaThreadGridTest);
TEST_COMPOSITE_ADD_UNIT(CudaComplexTest);
TEST_COMPOSITE_ADD_UNIT(CudaVecOpTest);
TEST_COMPOSITE_ADD_UNIT(CudaReduceTest);
TEST_COMPOSITE_END

#endif
