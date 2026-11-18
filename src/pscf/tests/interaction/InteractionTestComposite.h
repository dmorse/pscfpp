#ifndef PSCF_INTERACTION_TEST_COMPOSITE_H
#define PSCF_INTERACTION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "InteractionTest.h"
#include "AmbdInteractionTest.h"
#include "AmIteratorTmplTest.h"

TEST_COMPOSITE_BEGIN(InteractionTestComposite)
TEST_COMPOSITE_ADD_UNIT(InteractionTest);
TEST_COMPOSITE_ADD_UNIT(AmbdInteractionTest);
TEST_COMPOSITE_ADD_UNIT(AmIteratorTmplTest);
TEST_COMPOSITE_END

#endif
