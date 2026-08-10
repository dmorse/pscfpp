#ifndef RPG_TEST_SYSTEM_TEST_COMPOSITE_H
#define RPG_TEST_SYSTEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "SystemTest.h"
#include "AmIteratorBasisTest.h"
#include "AmIteratorGridTest.h"
#include "BeadTest.h"
#include "ThreadTest.h"

TEST_COMPOSITE_BEGIN(SystemTestComposite)
TEST_COMPOSITE_ADD_UNIT(SystemTest);
TEST_COMPOSITE_ADD_UNIT(AmIteratorBasisTest);
TEST_COMPOSITE_ADD_UNIT(AmIteratorGridTest);
TEST_COMPOSITE_ADD_UNIT(BeadTest)
TEST_COMPOSITE_ADD_UNIT(ThreadTest)
TEST_COMPOSITE_END

#endif
