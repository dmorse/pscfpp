#ifndef PSSP_BASIS_TEST_H
#define PSSP_BASIS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

#include <pssp/basis/Basis.h>
#include <pssp/basis/TWave.h>
#include <pssp/field/RField.h>
#include <pssp/field/RFieldDft.h>
#include <pssp/field/FFT.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/MeshIterator.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class BasisTest : public UnitTest 
{

public:

   void setUp()
   {
      setVerbose(2);
   }

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      Basis<3> sampleBasis;
      TEST_ASSERT(eq(sampleBasis.nWave(),0));
      TEST_ASSERT(eq(sampleBasis.nStar(),0));
      TEST_ASSERT(eq(sampleBasis.nBasis(),0));
   }

   void testMake1DBasis_I() {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<1> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_lamellar", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<1> d;
      d[0] = 8;
      Mesh<1> mesh(d);

      // Construct basis object using identity space group
      Basis<1> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);

      TEST_ASSERT(eq(basis.nWave(),8));
      TEST_ASSERT(eq(basis.nStar(),8));
      TEST_ASSERT(eq(basis.nBasis(),8));

      #if 0
      basis.outputWaves(std::cout);
      basis.outputStars(std::cout);
      #endif
   }

   void testMake1DBasis_Ib() {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<1> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_lamellar", in);
      in >> unitCell;

      // Make mesh object
      IntVec<1> d;
      in >> d;
      in.close();
      Mesh<1> mesh(d);

      // Read space group
      SpaceGroup<1> group;
      openInputFile("in/Group_1D_symmetric", in);
      in >> group;
      in.close();
      TEST_ASSERT(group.size() == 2);

      // Construct basis object using identity space group
      Basis<1> basis;
      basis.makeBasis(mesh, unitCell, group);

      TEST_ASSERT(eq(basis.nWave(), 8));
      TEST_ASSERT(eq(basis.nStar(), 5));
      TEST_ASSERT(eq(basis.nBasis(), 5));

      #if 0
      if (verbose() > 1) {
         basis.outputWaves(std::cout);
         basis.outputStars(std::cout);
      }
      #endif

   }

   void testMake2DBasis_square_33()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_square", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<2> d;
      d[0] = 3;
      d[1] = 3;
      Mesh<2> mesh(d);

      // Construct basis object using identity space group
      Basis<2> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);
   
      TEST_ASSERT(eq(basis.nWave(), 9));
      TEST_ASSERT(eq(basis.nStar(),9));

      #if 0
      basis.outputWaves(std::cout);   
      basis.outputStars(std::cout);   
      #endif

   }

   void testMake2DBasis_square_44()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_square", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<2> d;
      d[0] = 4;
      d[1] = 4;
      Mesh<2> mesh(d);

      // Read space group
      SpaceGroup<2> group;
      openInputFile("in/Group_2D_CenteredSquare", in);
      in >> group;

      // Make basis
      Basis<2> basis;
      basis.makeBasis(mesh, unitCell, group);
     
      TEST_ASSERT(eq(basis.nWave(), 16));
      TEST_ASSERT(eq(basis.nStar(), 6));
      TEST_ASSERT(eq(basis.nBasis(), 4));
      TEST_ASSERT(basis.isValid());

      #if 0
      if (verbose() > 1) {
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
      }
      #endif

   }

   void testMake2DBasis_hex()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Read UnitCell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_hexagonal", in);
      in >> unitCell;
      in.close();

      // Make Mesh object
      IntVec<2> d;
      d[0] = 24;
      d[1] = 24;
      Mesh<2> mesh(d);

      // Read space group
      SpaceGroup<2> group;
      openInputFile("in/p_6_m_m", in);
      in >> group;
      in.close();

      // Make basis
      Basis<2> basis;
      basis.makeBasis(mesh, unitCell, group);
     
      TEST_ASSERT(basis.isValid());
      //TEST_ASSERT(eq(basis.nWave(), 16));
      //TEST_ASSERT(eq(basis.nStar(), 6));
      //TEST_ASSERT(eq(basis.nBasis(), 4));

      #if 0
      if (verbose() > 1) {
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
      }
      #endif

   }

   void testMake3DBasis_I()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_cubic", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<3> d;
      d[0] = 3;
      d[1] = 3;
      d[2] = 3;
      Mesh<3> mesh(d);

      // Construct basis object
      Basis<3> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(basis.nWave(), 27));
      TEST_ASSERT(eq(basis.nStar(), 27));
      TEST_ASSERT(basis.isValid());
   }

   void testMake3DBasis_I_m_3b_m()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_cubic", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<3> d;
      d[0] = 8;
      d[1] = 8;
      d[2] = 8;
      Mesh<3> mesh(d);

      // Read group
      SpaceGroup<3> group;
      openInputFile("in/I_m_-3_m", in);
      in >> group;
      in.close();

      // Construct basis object
      Basis<3> basis;
      basis.makeBasis(mesh, unitCell, group);
      
      TEST_ASSERT(basis.isValid());
      //TEST_ASSERT(eq(basis.nWave(), 512));

      #if 0
      if (verbose() > 1) {
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
      }
      #endif

   }

   void testMake3DBasis_I_a_3b_d()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell_cubic", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<3> d;
      d[0] = 8;
      d[1] = 8;
      d[2] = 8;
      Mesh<3> mesh(d);

      // Read group
      SpaceGroup<3> group;
      openInputFile("in/I_a_-3_d", in);
      in >> group;
      in.close();

      // Construct basis object
      Basis<3> basis;
      basis.makeBasis(mesh, unitCell, group);
      
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));

      #if 0
      if (verbose() > 1) {
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
      }
      #endif

   }

};

TEST_BEGIN(BasisTest)
TEST_ADD(BasisTest, testConstructor)
TEST_ADD(BasisTest, testMake1DBasis_I)
TEST_ADD(BasisTest, testMake1DBasis_Ib)
TEST_ADD(BasisTest, testMake2DBasis_square_33)
TEST_ADD(BasisTest, testMake2DBasis_square_44)
TEST_ADD(BasisTest, testMake2DBasis_hex)
TEST_ADD(BasisTest, testMake3DBasis_I)
TEST_ADD(BasisTest, testMake3DBasis_I_m_3b_m)
TEST_ADD(BasisTest, testMake3DBasis_I_a_3b_d) 
TEST_END(BasisTest)

#endif
