#ifndef PSCF_TEST_CRYSTAL_TEST_COMPOSITE_H
#define PSCF_TEST_CRYSTAL_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "UnitCellTest.h"
#include "SpaceSymmetryTest.h"

TEST_COMPOSITE_BEGIN(CrystalTestComposite)
TEST_COMPOSITE_ADD_UNIT(UnitCellTest);
TEST_COMPOSITE_ADD_UNIT(SpaceSymmetryTest);
TEST_COMPOSITE_END

#endif
