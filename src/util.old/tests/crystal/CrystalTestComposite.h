#ifndef CRYSTAL_TEST_COMPOSITE_H
#define CRYSTAL_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "PointSymmetryTest.h"
#include "PointGroupTest.h"

TEST_COMPOSITE_BEGIN(CrystalTestComposite)
TEST_COMPOSITE_ADD_UNIT(PointSymmetryTest);
TEST_COMPOSITE_ADD_UNIT(PointGroupTest);
TEST_COMPOSITE_END

#endif
