#ifndef RPC_TEST_ITERATOR_TEST_COMPOSITE_H
#define RPC_TEST_ITERATOR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "MaskGenFilmTest.h"
#include "ExtGenFilmTest.h"
#include "ImposedFieldsGeneratorTest.h"

TEST_COMPOSITE_BEGIN(IteratorTestComposite)
TEST_COMPOSITE_ADD_UNIT(MaskGenFilmTest)
TEST_COMPOSITE_ADD_UNIT(ExtGenFilmTest)
TEST_COMPOSITE_ADD_UNIT(ImposedFieldsGeneratorTest)
TEST_COMPOSITE_END

#endif 
