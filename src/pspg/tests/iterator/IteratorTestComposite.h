#ifndef PSPG_TEST_ITERATOR_TEST_COMPOSITE_H
#define PSPG_TEST_ITERATOR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "IteratorMediatorCUDATest.h"
//#include "AmStrategyCUDATest.h"

TEST_COMPOSITE_BEGIN(IteratorTestComposite)
TEST_COMPOSITE_ADD_UNIT(IteratorMediatorCUDATest);
TEST_COMPOSITE_END

#endif
