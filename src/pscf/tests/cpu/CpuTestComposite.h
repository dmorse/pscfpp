#ifndef PSCF_CPU_TEST_COMPOSITE_H
#define PSCF_CPU_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CpuComplexTest.h"
#include "CpuVecOpTest.h"
#include "CpuVecRandomTest.h"
#include "CpuFftwDArrayTest.h"
#include "CpuFftwDRArrayTest.h"
#include "CpuDeviceArrayTest.h"
#include "CpuHostArrayTest.h"
#include "CpuConstHostArrayTest.h"

TEST_COMPOSITE_BEGIN(CpuTestComposite)
TEST_COMPOSITE_ADD_UNIT(CpuComplexTest);
TEST_COMPOSITE_ADD_UNIT(CpuVecOpTest);
TEST_COMPOSITE_ADD_UNIT(CpuVecRandomTest);
TEST_COMPOSITE_ADD_UNIT(CpuFftwDArrayTest);
TEST_COMPOSITE_ADD_UNIT(CpuFftwDRArrayTest);
TEST_COMPOSITE_ADD_UNIT(CpuDeviceArrayTest);
TEST_COMPOSITE_ADD_UNIT(CpuHostArrayTest);
TEST_COMPOSITE_ADD_UNIT(CpuConstHostArrayTest);
TEST_COMPOSITE_END

#endif
