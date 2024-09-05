#ifndef RPC_SIMULATE_TEST_COMPOSITE_H
#define RPC_SIMULATE_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "SimulatorTest.h"
#include "CompressorTest.h"
#include "McSimulatorTest.h"
#include "BdSimulatorTest.h"

TEST_COMPOSITE_BEGIN(SimulateTestComposite)
TEST_COMPOSITE_ADD_UNIT(SimulatorTest)
TEST_COMPOSITE_ADD_UNIT(CompressorTest)
TEST_COMPOSITE_ADD_UNIT(McSimulatorTest)
TEST_COMPOSITE_ADD_UNIT(BdSimulatorTest)
TEST_COMPOSITE_END

#endif 
