#ifndef PSCF_FIELD_CUDA_TEST_COMPOSITE_H
#define PSCF_FIELD_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaFieldTest.h"
#include "CudaFieldComparisonTest.h"
#include "CudaFftTest.h"
#include "CudaWaveListTest.h"

TEST_COMPOSITE_BEGIN(FieldCudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaFieldTest);
TEST_COMPOSITE_ADD_UNIT(CudaFieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(CudaFftTest);
TEST_COMPOSITE_ADD_UNIT(CudaWaveListTest);
TEST_COMPOSITE_END

#endif
