#ifndef PRDC_TEST_CRYSTAL_TEST_COMPOSITE_H
#define PRDC_TEST_CRYSTAL_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "UnitCellTest.h"
#include "SpaceSymmetryTest.h"
#include "SpaceGroupTest.h"
#include "BasisTest.h"

TEST_COMPOSITE_BEGIN(CrystalTestComposite)
TEST_COMPOSITE_ADD_UNIT(UnitCellTest);
TEST_COMPOSITE_ADD_UNIT(SpaceSymmetryTest);
TEST_COMPOSITE_ADD_UNIT(SpaceGroupTest);
TEST_COMPOSITE_ADD_UNIT(BasisTest);
TEST_COMPOSITE_END

#endif
