#ifndef PSCF_SOLVERS_TEST_COMPOSITE_H
#define PSCF_SOLVERS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "BlockStubTest.h"
#include "PolymerStubTest.h"
#include "MixtureStubTest.h"

TEST_COMPOSITE_BEGIN(SolversTestComposite)
TEST_COMPOSITE_ADD_UNIT(BlockStubTest);
TEST_COMPOSITE_ADD_UNIT(PolymerStubTest);
TEST_COMPOSITE_ADD_UNIT(MixtureStubTest);
TEST_COMPOSITE_END

#endif
