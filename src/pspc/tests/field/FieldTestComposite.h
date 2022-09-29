#ifndef PSPC_FIELD_TEST_COMPOSITE_H
#define PSPC_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldTest.h"
#include "RFieldTest.h"
#include "RFieldDftTest.h"
#include "FftTest.h"
#include "FieldComparisonTest.h"
#include "DomainTest.h"
#include "FieldIoTest.h"
#include "FieldContainerTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldTest);
TEST_COMPOSITE_ADD_UNIT(RFieldDftTest);
TEST_COMPOSITE_ADD_UNIT(FftTest);
TEST_COMPOSITE_ADD_UNIT(FieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoTest);
TEST_COMPOSITE_ADD_UNIT(FieldContainerTest);
TEST_COMPOSITE_END

#endif
