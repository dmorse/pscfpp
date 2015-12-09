#ifndef FD1D_TEST_COMPOSITE_H
#define FD1D_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "TridiagonalSolverTest.h"
#include "PropagatorTest.h"

TEST_COMPOSITE_BEGIN(Fd1dTestComposite)
TEST_COMPOSITE_ADD_UNIT(TridiagonalSolverTest);
TEST_COMPOSITE_ADD_UNIT(PropagatorTest);
TEST_COMPOSITE_END

#endif
