#ifndef RPC_TEST_SYSTEM_TEST_COMPOSITE_H
#define RPC_TEST_SYSTEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "SystemTest.h"
#include "BeadTest.h"
#include "ThreadTest.h"

TEST_COMPOSITE_BEGIN(SystemTestComposite)
TEST_COMPOSITE_ADD_UNIT(SystemTest)
TEST_COMPOSITE_ADD_UNIT(BeadTest)
TEST_COMPOSITE_ADD_UNIT(ThreadTest)
TEST_COMPOSITE_END

#endif 
