#ifndef MPI_TEST_COMPOSITE_H
#define MPI_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MpiSendRecvTest.h"
#include "MpiFileIoTest.h"
#include "MpiLoaderTest.h"
//#include "MpiLoggerTest.h"

using namespace Util;

TEST_COMPOSITE_BEGIN(MpiTestComposite)

TEST_COMPOSITE_ADD_UNIT(MpiSendRecvTest)
TEST_COMPOSITE_ADD_UNIT(MpiFileIoTest)
TEST_COMPOSITE_ADD_UNIT(MpiLoaderTest)

TEST_COMPOSITE_END

#endif
