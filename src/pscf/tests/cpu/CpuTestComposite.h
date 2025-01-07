#ifndef PSCF_TEST_CPU_TEST_COMPOSITE_H
#define PSCF_TEST_CPU_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CpuComplexTest.h"

TEST_COMPOSITE_BEGIN(CpuTestComposite)
TEST_COMPOSITE_ADD_UNIT(CpuComplexTest);
TEST_COMPOSITE_END

#endif
