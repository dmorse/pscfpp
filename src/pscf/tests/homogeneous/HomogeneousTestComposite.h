#ifndef PSCF_HOMOGENEOUS_TEST_COMPOSITE_H
#define PSCF_HOMOGENEOUS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "GroupTest.h"

TEST_COMPOSITE_BEGIN(HomogeneousTestComposite)
TEST_COMPOSITE_ADD_UNIT(GroupTest);
TEST_COMPOSITE_END

#endif
