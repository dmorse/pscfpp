#ifndef PSSP_FIELD_TEST_COMPOSITE_H
#define PSSP_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldTest.h"
#include "RFieldTest.h"
#include "RFieldDftTest.h"
#include "FftTest.h"
//#include "FieldUtilTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldDftTest);
TEST_COMPOSITE_ADD_UNIT(FftTest);
//TEST_COMPOSITE_ADD_UNIT(FieldUtilTest);
TEST_COMPOSITE_END

#endif
