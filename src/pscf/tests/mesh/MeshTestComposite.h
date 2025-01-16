#ifndef PSCF_TEST_MESH_TEST_COMPOSITE_H
#define PSCF_TEST_MESH_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MeshTest.h"
#include "MeshIteratorTest.h"
#include "MeshIteratorFortranTest.h"

TEST_COMPOSITE_BEGIN(MeshTestComposite)
TEST_COMPOSITE_ADD_UNIT(MeshTest);
TEST_COMPOSITE_ADD_UNIT(MeshIteratorTest);
TEST_COMPOSITE_ADD_UNIT(MeshIteratorFortranTest);
TEST_COMPOSITE_END

#endif
