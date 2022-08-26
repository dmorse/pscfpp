#ifndef PSPG_TEST_SYSTEM_TEST_COMPOSITE_H
#define PSPG_TEST_SYSTEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "SystemTest.h"
#include "AmIteratorTest.h"

TEST_COMPOSITE_BEGIN(SystemTestComposite)
TEST_COMPOSITE_ADD_UNIT(SystemTest);
TEST_COMPOSITE_ADD_UNIT(AmIteratorTest);
TEST_COMPOSITE_END

#endif
