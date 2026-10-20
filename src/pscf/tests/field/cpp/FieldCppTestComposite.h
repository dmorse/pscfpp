#ifndef PSCF_FIELD_CPP_TEST_COMPOSITE_H
#define PSCF_FIELD_CPP_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CppRFieldTest.h"
#include "CppRFieldDftTest.h"
#include "CppCFieldTest.h"
#include "CppFieldComparisonTest.h"
#include "CppFftTest.h"
#include "CppWaveListTest.h"

TEST_COMPOSITE_BEGIN(FieldCppTestComposite)
TEST_COMPOSITE_ADD_UNIT(CppRFieldTest);
TEST_COMPOSITE_ADD_UNIT(CppRFieldDftTest);
TEST_COMPOSITE_ADD_UNIT(CppCFieldTest);
TEST_COMPOSITE_ADD_UNIT(CppFieldComparisonTest);
TEST_COMPOSITE_ADD_UNIT(CppFftTest);
TEST_COMPOSITE_ADD_UNIT(CppWaveListTest);
TEST_COMPOSITE_END

#endif
