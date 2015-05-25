#ifndef CONTAINERS_TEST_COMPOSITE_H
#define CONTAINERS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "DArrayTest.h"
#include "FArrayTest.h"
#include "RArrayTest.h"

#include "DSArrayTest.h"
#include "FSArrayTest.h"
#include "GArrayTest.h"

#include "DPArrayTest.h"
#include "FPArrayTest.h"
#include "GPArrayTest.h"
#include "ArraySetTest.h"
#include "ArrayStackTest.h"
#include "RingBufferTest.h"
#include "SSetTest.h"

#include "DMatrixTest.h"
#include "FMatrixTest.h"
#include "DRaggedMatrixTest.h"

#include "GridArrayTest.h"

#include "ListTest.h"

TEST_COMPOSITE_BEGIN(ContainersTestComposite)

TEST_COMPOSITE_ADD_UNIT(DArrayTest)
TEST_COMPOSITE_ADD_UNIT(FArrayTest)
TEST_COMPOSITE_ADD_UNIT(RArrayTest)
TEST_COMPOSITE_ADD_UNIT(DSArrayTest)
TEST_COMPOSITE_ADD_UNIT(FSArrayTest)
TEST_COMPOSITE_ADD_UNIT(GArrayTest)

TEST_COMPOSITE_ADD_UNIT(FPArrayTest)
TEST_COMPOSITE_ADD_UNIT(DPArrayTest)
TEST_COMPOSITE_ADD_UNIT(GPArrayTest)
TEST_COMPOSITE_ADD_UNIT(ArraySetTest)
TEST_COMPOSITE_ADD_UNIT(ArrayStackTest)
TEST_COMPOSITE_ADD_UNIT(RingBufferTest)
TEST_COMPOSITE_ADD_UNIT(SSetTest)

TEST_COMPOSITE_ADD_UNIT(DMatrixTest)
TEST_COMPOSITE_ADD_UNIT(FMatrixTest)
TEST_COMPOSITE_ADD_UNIT(DRaggedMatrixTest)

TEST_COMPOSITE_ADD_UNIT(GridArrayTest)

TEST_COMPOSITE_ADD_UNIT(ListTest)

TEST_COMPOSITE_END

#endif
