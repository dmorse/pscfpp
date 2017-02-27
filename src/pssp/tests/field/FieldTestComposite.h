#ifndef PSSP_TEST_FIELD_TEST_COMPOSITE_H
#define PSSP_TEST_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FieldTest.h"
#include "RMeshFieldTest.h"
#include "KMeshFieldTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(FieldTest);
TEST_COMPOSITE_ADD_UNIT(RMeshFieldTest);
TEST_COMPOSITE_ADD_UNIT(KMeshFieldTest);
TEST_COMPOSITE_END

#endif
