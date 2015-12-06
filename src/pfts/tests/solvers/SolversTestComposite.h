#ifndef SOLVERS_TEST_COMPOSITE_H
#define SOLVERS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MonomerTest.h"

TEST_COMPOSITE_BEGIN(SolversTestComposite)
TEST_COMPOSITE_ADD_UNIT(MonomerTest);
TEST_COMPOSITE_END

#endif
