#ifndef PSCF_CPP_TEST_COMPOSITE_H
#define PSCF_CPP_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CppComplexTest.h"
#include "CppVecOpTest.h"
#include "CppVecRandomTest.h"
#include "CppFftwDArrayTest.h"
#include "CppFftwDRArrayTest.h"
#include "CppDeviceArrayTest.h"
#include "CppHostArrayTest.h"
#include "CppConstHostArrayTest.h"

TEST_COMPOSITE_BEGIN(CppTestComposite)
TEST_COMPOSITE_ADD_UNIT(CppComplexTest);
TEST_COMPOSITE_ADD_UNIT(CppVecOpTest);
TEST_COMPOSITE_ADD_UNIT(CppVecRandomTest);
TEST_COMPOSITE_ADD_UNIT(CppFftwDArrayTest);
TEST_COMPOSITE_ADD_UNIT(CppFftwDRArrayTest);
TEST_COMPOSITE_ADD_UNIT(CppDeviceArrayTest);
TEST_COMPOSITE_ADD_UNIT(CppHostArrayTest);
TEST_COMPOSITE_ADD_UNIT(CppConstHostArrayTest);
TEST_COMPOSITE_END

#endif
