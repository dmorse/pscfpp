#ifndef RPG_ENVIRONMENT_TEST_COMPOSITE_H
#define RPG_ENVIRONMENT_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "FilmFieldGenMaskTest.h"
#include "FilmFieldGenExtTest.h"
#include "FilmEnvironmentTest.h"

TEST_COMPOSITE_BEGIN(EnvironmentTestComposite)
TEST_COMPOSITE_ADD_UNIT(FilmFieldGenMaskTest);
TEST_COMPOSITE_ADD_UNIT(FilmFieldGenExtTest);
TEST_COMPOSITE_ADD_UNIT(FilmEnvironmentTest);
TEST_COMPOSITE_END

#endif
