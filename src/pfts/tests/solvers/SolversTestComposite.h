#ifndef SOLVERS_TEST_COMPOSITE_H
#define SOLVERS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "PolymerStubTest.h"
#include "SystemStubTest.h"

TEST_COMPOSITE_BEGIN(SolversTestComposite)
TEST_COMPOSITE_ADD_UNIT(PolymerStubTest);
TEST_COMPOSITE_ADD_UNIT(SystemStubTest);
TEST_COMPOSITE_END

#endif
