#ifndef ACCUMULATOR_TEST_COMPOSITE_H
#define ACCUMULATOR_TEST_COMPOSITE_H

#include "AverageTest.h"
#include "AutoCorrTest.h"
#include "AutoCorrArrayTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(AccumulatorTestComposite)
TEST_COMPOSITE_ADD_UNIT(AverageTest)
TEST_COMPOSITE_ADD_UNIT(AutoCorrTest)
TEST_COMPOSITE_ADD_UNIT(AutoCorrArrayTest)
TEST_COMPOSITE_END

#endif
