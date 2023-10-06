#ifndef PSPC_FIELD_TEST_COMPOSITE_H
#define PSPC_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldIoTest.h"
#include "DomainTest.h"
#include "WFieldContainerTest.h"
#include "CFieldContainerTest.h"
#include "MaskTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldIoTest);
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(WFieldContainerTest);
TEST_COMPOSITE_ADD_UNIT(CFieldContainerTest);
TEST_COMPOSITE_ADD_UNIT(MaskTest);
TEST_COMPOSITE_END

#endif
