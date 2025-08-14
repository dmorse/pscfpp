#ifndef PSCF_FLORY_HUGGINS_TEST_COMPOSITE_H
#define PSCF_FLORY_HUGGINS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ClumpTest.h"
#include "MoleculeTest.h"
#include "MixtureTest.h"

TEST_COMPOSITE_BEGIN(FloryHugginsTestComposite)
TEST_COMPOSITE_ADD_UNIT(ClumpTest);
TEST_COMPOSITE_ADD_UNIT(MoleculeTest);
TEST_COMPOSITE_ADD_UNIT(MixtureTest);
TEST_COMPOSITE_END

#endif
