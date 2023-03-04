#ifndef CYLN_FIELD_TEST_COMPOSITE_H
#define CYLN_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FftTest.h"
#include "FieldTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FftTest);
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_END

#endif
