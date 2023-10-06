#ifndef PRDC_CPU_TEST_COMPOSITE_H
#define PRDC_CPU_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldTest.h"
#include "RFieldTest.h"
#include "RFieldDftTest.h"
#include "FieldComparisonTest.h"
#include "FftTest.h"

TEST_COMPOSITE_BEGIN(CpuTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldDftTest);
TEST_COMPOSITE_ADD_UNIT(FieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(FftTest);
TEST_COMPOSITE_END

#endif
