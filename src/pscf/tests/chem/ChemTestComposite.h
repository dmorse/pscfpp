#ifndef PSCF_TESTS_CHEM_TEST_COMPOSITE_H
#define PSCF_TESTS_CHEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MonomerTest.h"
#include "VertexTest.h"
#include "BlockDescriptorTest.h"
#include "PolymerTypeTest.h"
#include "PolymerModelTest.h"

TEST_COMPOSITE_BEGIN(ChemTestComposite)
TEST_COMPOSITE_ADD_UNIT(MonomerTest);
TEST_COMPOSITE_ADD_UNIT(VertexTest);
TEST_COMPOSITE_ADD_UNIT(BlockDescriptorTest);
TEST_COMPOSITE_ADD_UNIT(PolymerTypeTest);
TEST_COMPOSITE_ADD_UNIT(PolymerModelTest);
TEST_COMPOSITE_END

#endif
