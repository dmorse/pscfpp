#ifndef MISC_TEST_COMPOSITE_H
#define MISC_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "SetableTest.h"
#include "BitTest.h"
#include "ioUtilTest.h"
#include "XmlTest.h"

TEST_COMPOSITE_BEGIN(MiscTestComposite)
TEST_COMPOSITE_ADD_UNIT(SetableTest);
TEST_COMPOSITE_ADD_UNIT(BitTest);
TEST_COMPOSITE_ADD_UNIT(ioUtilTest);
TEST_COMPOSITE_ADD_UNIT(XmlTest);
TEST_COMPOSITE_END

#endif
