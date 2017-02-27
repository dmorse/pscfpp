#ifndef PSSP_TEST_FIELD_TEST_COMPOSITE_H
#define PSSP_TEST_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldTest.h"
#include "RFieldTest.h"
#include "RFieldDFTTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldDFTTest);
TEST_COMPOSITE_END

#endif
