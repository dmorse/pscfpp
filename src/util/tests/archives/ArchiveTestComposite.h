#ifndef ARCHIVE_TEST_COMPOSITE_H
#define ARCHIVE_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MemoryArchiveTest.h"
#include "BinaryFileArchiveTest.h"
#include "TextFileArchiveTest.h"
#include "XdrFileArchiveTest.h"

TEST_COMPOSITE_BEGIN(ArchiveTestComposite)
TEST_COMPOSITE_ADD_UNIT(MemoryArchiveTest);
TEST_COMPOSITE_ADD_UNIT(BinaryFileArchiveTest);
TEST_COMPOSITE_ADD_UNIT(TextFileArchiveTest);
TEST_COMPOSITE_ADD_UNIT(XdrFileArchiveTest);
TEST_COMPOSITE_END

#endif
