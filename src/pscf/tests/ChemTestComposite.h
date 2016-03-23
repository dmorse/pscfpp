#ifndef PSCF_CHEM_TEST_COMPOSITE_H
#define PSCF_CHEM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MonomerTest.h"
#include "VertexTest.h"
#include "BlockDescriptorTest.h"
#include "BlockStubTest.h"
#include "PolymerStubTest.h"
#include "MixtureStubTest.h"
#include "ChiInteractionTest.h"
#include "TridiagonalSolverTest.h"
#ifdef PSCF_GSL
#include "LuSolverTest.h"
#endif
#include "IntVecTest.h"

TEST_COMPOSITE_BEGIN(ChemTestComposite)
TEST_COMPOSITE_ADD_UNIT(MonomerTest);
TEST_COMPOSITE_ADD_UNIT(VertexTest);
TEST_COMPOSITE_ADD_UNIT(BlockDescriptorTest);
TEST_COMPOSITE_ADD_UNIT(BlockStubTest);
TEST_COMPOSITE_ADD_UNIT(PolymerStubTest);
TEST_COMPOSITE_ADD_UNIT(MixtureStubTest);
TEST_COMPOSITE_ADD_UNIT(ChiInteractionTest);
TEST_COMPOSITE_ADD_UNIT(TridiagonalSolverTest);
#ifdef PSCF_GSL
TEST_COMPOSITE_ADD_UNIT(LuSolverTest);
#endif
TEST_COMPOSITE_ADD_UNIT(IntVecTest);
TEST_COMPOSITE_END

#endif
