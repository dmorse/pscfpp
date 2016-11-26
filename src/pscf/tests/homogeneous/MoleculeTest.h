#ifndef PSCF_HOMOGENEOUS_MOLECULE_TEST_H
#define PSCF_HOMOGENEOUS_MOLECULE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/homogeneous/Molecule.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class MoleculeTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Homogeneous::Molecule molecule;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Homogeneous::Molecule molecule;
      std::ifstream in;
      openInputFile("in/Molecule", in);

      molecule.readParam(in);
      TEST_ASSERT(molecule.nClump() == 2);
      TEST_ASSERT(molecule.clump(0).monomerId() == 0);
      TEST_ASSERT(eq(molecule.clump(0).size(), 2.0));
      TEST_ASSERT(molecule.clump(1).monomerId() == 1);
      TEST_ASSERT(eq(molecule.clump(1).size(), 3.0));
      TEST_ASSERT(eq(molecule.size(), 5.0));
      molecule.writeParam(std::cout) ;
   }

   void testSetters()
   {
      printMethod(TEST_FUNC);
      Homogeneous::Molecule molecule;

      molecule.setNClump(2);
      molecule.clump(0).setMonomerId(0);
      molecule.clump(0).setSize(2.0);
      molecule.clump(1).setMonomerId(1);
      molecule.clump(1).setSize(3.0);
      molecule.computeSize();

      TEST_ASSERT(molecule.nClump() == 2);
      TEST_ASSERT(molecule.clump(0).monomerId() == 0);
      TEST_ASSERT(eq(molecule.clump(0).size(), 2.0));
      TEST_ASSERT(molecule.clump(1).monomerId() == 1);
      TEST_ASSERT(eq(molecule.clump(1).size(), 3.0));
      TEST_ASSERT(eq(molecule.size(), 5.0));
   } 

};

TEST_BEGIN(MoleculeTest)
TEST_ADD(MoleculeTest, testConstructor)
TEST_ADD(MoleculeTest, testReadWrite)
TEST_ADD(MoleculeTest, testSetters)
TEST_END(MoleculeTest)

#endif
