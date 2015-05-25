#ifndef SPACE_TEST_COMPOSITE_H
#define SPACE_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "VectorTest.h"
#include "IntVectorTest.h"
#include "TensorTest.h"
#include "GridTest.h"

TEST_COMPOSITE_BEGIN(SpaceTestComposite)
TEST_COMPOSITE_ADD_UNIT(VectorTest);
TEST_COMPOSITE_ADD_UNIT(IntVectorTest);
TEST_COMPOSITE_ADD_UNIT(TensorTest);
TEST_COMPOSITE_ADD_UNIT(GridTest);
TEST_COMPOSITE_END

#endif
