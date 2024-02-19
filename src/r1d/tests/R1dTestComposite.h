#ifndef R1D_TEST_COMPOSITE_H
#define R1D_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "DomainTest.h"
#include "PropagatorTest.h"
#include "MixtureTest.h"
#include "SystemTest.h"

TEST_COMPOSITE_BEGIN(R1dTestComposite)
TEST_COMPOSITE_ADD_UNIT(DomainTest);
TEST_COMPOSITE_ADD_UNIT(PropagatorTest);
TEST_COMPOSITE_ADD_UNIT(MixtureTest);
TEST_COMPOSITE_ADD_UNIT(SystemTest);
TEST_COMPOSITE_END

#endif
