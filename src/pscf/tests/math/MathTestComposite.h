#ifndef PSCF_TEST_MATH_TEST_COMPOSITE_H
#define PSCF_TEST_MATH_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "IntVecTest.h"
#include "RealVecTest.h"
#include "TridiagonalSolverTest.h"
#include "LuSolverTest.h"

TEST_COMPOSITE_BEGIN(MathTestComposite)
TEST_COMPOSITE_ADD_UNIT(IntVecTest);
TEST_COMPOSITE_ADD_UNIT(RealVecTest);
TEST_COMPOSITE_ADD_UNIT(TridiagonalSolverTest);
TEST_COMPOSITE_ADD_UNIT(LuSolverTest);
TEST_COMPOSITE_END

#endif
