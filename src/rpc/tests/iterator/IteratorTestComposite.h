#ifndef RPC_TEST_ITERATOR_TEST_COMPOSITE_H
#define RPC_TEST_ITERATOR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "FilmIteratorTest.h"

TEST_COMPOSITE_BEGIN(IteratorTestComposite)
TEST_COMPOSITE_ADD_UNIT(FilmIteratorTest)
TEST_COMPOSITE_END

#endif 
