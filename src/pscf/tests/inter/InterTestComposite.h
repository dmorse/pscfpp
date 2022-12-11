#ifndef PSCF_INTER_TEST_COMPOSITE_H
#define PSCF_INTER_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "InteractionTest.h"
#include "AmbdInteractionTest.h"

TEST_COMPOSITE_BEGIN(InterTestComposite)
TEST_COMPOSITE_ADD_UNIT(InteractionTest);
TEST_COMPOSITE_ADD_UNIT(AmbdInteractionTest);
TEST_COMPOSITE_END

#endif
