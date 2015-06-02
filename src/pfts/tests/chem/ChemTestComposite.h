#ifndef CHEM_TEST_COMPOSITE_H
#define CHEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MonomerTest.h"
#include "BlockTest.h"
#include "VertexTest.h"
#include "PolymerTest.h"

TEST_COMPOSITE_BEGIN(ChemTestComposite)
TEST_COMPOSITE_ADD_UNIT(MonomerTest);
TEST_COMPOSITE_ADD_UNIT(BlockTest);
TEST_COMPOSITE_ADD_UNIT(VertexTest);
TEST_COMPOSITE_ADD_UNIT(PolymerTest);
TEST_COMPOSITE_END

#endif
