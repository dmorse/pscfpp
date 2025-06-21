#ifndef RPG_FIELD_TEST_COMPOSITE_H
#define RPG_FIELD_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "DomainTest.h"
#include "FieldIoTest.h"
#include "MaskTest.h"
#include "FilmFieldGenMaskTest.h"
#include "FilmFieldGenExtTest.h"
#include "FilmEnvironmentTest.h"

TEST_COMPOSITE_BEGIN(FieldTestComposite)
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(FieldIoTest);
TEST_COMPOSITE_ADD_UNIT(MaskTest);
TEST_COMPOSITE_ADD_UNIT(FilmFieldGenMaskTest);
TEST_COMPOSITE_ADD_UNIT(FilmFieldGenExtTest);
TEST_COMPOSITE_ADD_UNIT(FilmEnvironmentTest);
TEST_COMPOSITE_END

#endif
