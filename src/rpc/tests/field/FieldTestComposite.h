#ifndef RPC_FIELD_TEST_COMPOSITE_H
#define RPC_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldIoTest.h"
#include "DomainTest.h"
#include "WFieldsTest.h"
#include "CFieldsTest.h"
#include "MaskTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldIoTest);
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(WFieldsTest);
TEST_COMPOSITE_ADD_UNIT(CFieldsTest);
TEST_COMPOSITE_ADD_UNIT(MaskTest);
TEST_COMPOSITE_END

#endif
