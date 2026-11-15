#ifndef PSSP_BASIS_TEST_H
#define PSSP_BASIS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/BWave.h>
#include <pscf/crystal/groupFile.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>
#include <util/misc/Exception.h>

#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

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
      openInputFile("in/Lamellar", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<1> d;
      d[0] = 8;
      Mesh<1> mesh(d);

      // Construct basis object using identity space group
      Basis<1> basis;
      TEST_ASSERT(!basis.isInitialized());
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);

      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(eq(basis.nWave(), 8));
      TEST_ASSERT(eq(basis.nStar(), 8));
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
      openInputFile("in/Lamellar", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<1> d;
      d[0] = 8;
      Mesh<1> mesh(d);

      // Read space group
      SpaceGroup<1> group;
      openInputFile("in/Group_1D_Symmetric", in);
      in >> group;
      in.close();
      TEST_ASSERT(group.size() == 2);

      // Construct basis object using identity space group
      Basis<1> basis;
      basis.makeBasis(mesh, unitCell, group);

      TEST_ASSERT(basis.isInitialized());
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

   /*
   * Write the wave and star lists of a basis on a 1D mesh of 20001
   * points, with wave indices from -10000 to 10000, and check that
   * every integer column is a separate word with the right value.
   */
   void testOutputWavesStars_1D()
   {
      printMethod(TEST_FUNC);

      UnitCell<1> unitCell;
      std::ifstream in;
      openInputFile("in/Lamellar", in);
      in >> unitCell;
      in.close();
      IntVec<1> d;
      d[0] = 20001;
      Mesh<1> mesh(d);
      Basis<1> basis;
      basis.makeBasis(mesh, unitCell, "I");
      TEST_ASSERT(eq(basis.nStar(), 20001));

      // Split a line into words separated by white space
      auto split = [](std::string const & line) -> std::vector<std::string> {
         std::istringstream words(line);
         std::vector<std::string> list;
         std::string word;
         while (words >> word) {
            list.push_back(word);
         }
         return list;
      };

      // Waves: id, id, index, star id, real and imaginary coefficient
      std::stringstream waves;
      basis.outputWaves(waves);
      std::string line;
      std::getline(waves, line);
      std::getline(waves, line);
      for (int i = 0; i < basis.nWave(); ++i) {
         std::getline(waves, line);
         std::vector<std::string> words = split(line);
         TEST_ASSERT(eq((int)words.size(), 6));
         TEST_ASSERT(eq(std::stoi(words[1]), i));
         TEST_ASSERT(eq(std::stoi(words[2]), basis.wave(i).indicesMin[0]));
         TEST_ASSERT(eq(std::stoi(words[3]), basis.wave(i).starId));
      }

      // Stars: basis id, star id, size, begin, end, invertFlag, index
      std::stringstream stars;
      basis.outputStars(stars);
      std::getline(stars, line);
      std::getline(stars, line);
      for (int i = 0; i < basis.nStar(); ++i) {
         std::getline(stars, line);
         std::vector<std::string> words = split(line);
         TEST_ASSERT(eq((int)words.size(), 7));
         TEST_ASSERT(eq(std::stoi(words[1]), i));
         TEST_ASSERT(eq(std::stoi(words[6]), basis.star(i).waveMin[0]));
      }
   }

   void testMake2DBasis_square_33()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/Square", in);
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
   
      TEST_ASSERT(basis.isInitialized());
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
      openInputFile("in/Square", in);
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
     
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(eq(basis.nWave(), 16));
      TEST_ASSERT(eq(basis.nStar(), 6));
      TEST_ASSERT(eq(basis.nBasis(), 4));
      TEST_ASSERT(basis.isValid());

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/Square_44", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
      }

   }

   /*
   * Make a 2D basis. Return the number of basis functions, or -1 if an
   * Exception is thrown.
   */
   int nBasis2D(std::string const & cell, int n, std::string groupName)
   {
      UnitCell<2> unitCell;
      std::istringstream in(cell);
      in >> unitCell;
      IntVec<2> d;
      d[0] = n;
      d[1] = n;
      Mesh<2> mesh(d);
      Basis<2> basis;
      try {
         basis.makeBasis(mesh, unitCell, groupName);
      } catch (Exception& e) {
         return -1;
      }
      return basis.nBasis();
   }

   /*
   * Squared norms of waves are compared with a tolerance relative to
   * their magnitude, so that bases can be made for small cells and for
   * cell parameters rounded to the 11 significant digits of field
   * files.
   */
   void testMake2DBasis_normTolerance()
   {
      printMethod(TEST_FUNC);

      // Hexagonal cells of different sizes, p_6_m_m on a 24 x 24 mesh
      TEST_ASSERT(eq(nBasis2D("hexagonal 1.0", 24, "p_6_m_m"), 61));
      TEST_ASSERT(eq(nBasis2D("hexagonal 0.01", 24, "p_6_m_m"), 61));
      TEST_ASSERT(eq(nBasis2D("hexagonal 1.0E-4", 24, "p_6_m_m"), 61));

      // Rhombic cell with gamma = 2 pi/3, exact and rounded
      std::string exact = "rhombic 1.7 2.0943951023931953";
      std::string rounded = "rhombic 1.7 2.0943951024";
      TEST_ASSERT(eq(nBasis2D(exact, 16, "p_2"), 130));
      TEST_ASSERT(eq(nBasis2D(rounded, 16, "p_2"), 130));
      TEST_ASSERT(eq(nBasis2D(exact, 16, "p_6_m_m"), 30));
      TEST_ASSERT(eq(nBasis2D(rounded, 16, "p_6_m_m"), 30));
      TEST_ASSERT(eq(nBasis2D(exact, 13, "p_6_m_m"), 21));
      TEST_ASSERT(eq(nBasis2D(rounded, 13, "p_6_m_m"), 21));

      // A rhombic angle that does not have hexagonal symmetry
      std::string other = "rhombic 1.7 2.0944051024";
      TEST_ASSERT(eq(nBasis2D(other, 16, "p_6_m_m"), -1));
   }

   void testMake2DBasis_hex()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Read UnitCell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/Hexagonal2D", in);
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
     
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      //TEST_ASSERT(eq(basis.nWave(), 16));
      //TEST_ASSERT(eq(basis.nStar(), 6));
      //TEST_ASSERT(eq(basis.nBasis(), 4));

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/p_6_m_m", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
      }

   }

   void testMake3DBasis_I()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/Cubic", in);
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
      
      TEST_ASSERT(basis.isInitialized());
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
      openInputFile("in/Cubic", in);
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
      
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));

      // Check basisFunction
      int i, j;
      j = 0;
      for (i = 0; i < basis.nStar(); ++i) {
         if (!(basis.star(i).cancel)) {
            TEST_ASSERT(&basis.star(i) == &basis.basisFunction(j));
            ++j;
         }
      }
      TEST_ASSERT(j == basis.nBasis());

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/I_m_3b_m", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
      }

   }

   void testMake3DBasis_I_a_3b_d()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/Cubic", in);
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
      
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));

      // Compare star and basisFunction accessors
      int i, j;
      j = 0;
      for (i = 0; i < basis.nStar(); ++i) {
         if (!(basis.star(i).cancel)) {
            TEST_ASSERT(&basis.star(i) == &basis.basisFunction(j));
            ++j;
         }
      }
      TEST_ASSERT(j == basis.nBasis());

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/I_a_3b_d", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
      }

   }

   void testMake3DBasis_F_d_3b_m_2()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/Cubic", in);
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
      openInputFile("in/F_d_-3_m:2", in);
      // This is a centro-symmetric setting for the C15 Laves phase
      in >> group;
      in.close();


      // Construct basis object
      Basis<3> basis;
      basis.makeBasis(mesh, unitCell, group);
      
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));

      // Compare star and basisFunction accessors
      int i, j;
      j = 0;
      for (i = 0; i < basis.nStar(); ++i) {
         if (!(basis.star(i).cancel)) {
            TEST_ASSERT(&basis.star(i) == &basis.basisFunction(j));
            ++j;
         }
      }
      TEST_ASSERT(j == basis.nBasis());

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/F_d_-3_m:2", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
         out.close();
      }

   }

   void testMake3DBasis_F_d_3b_m_1()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/Cubic", in);
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
      openInputFile("in/F_d_-3_m:1", in);
      // This is a non-centro-symmetric setting for the C15 Laves phase
      in >> group;
      in.close();

      // Construct basis object
      Basis<3> basis;
      basis.makeBasis(mesh, unitCell, group);
  
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 512));

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/F_d_-3_m:1", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
         out.close();
      }

   }

   void testMake3DBasis_I_41_3_2()
   {
      printMethod(TEST_FUNC);
      // printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/Cubic", in);
      in >> unitCell;
      in.close();

      // Make mesh object
      IntVec<3> d;
      d[0] = 16;
      d[1] = 16;
      d[2] = 16;
      Mesh<3> mesh(d);

      // Read group
      SpaceGroup<3> group;
      openInputFile("in/I_41_3_2", in);
      // This is a non-centro-symmetric alternating gyroid group
      in >> group;
      in.close();

      // Construct basis object
      Basis<3> basis;
      basis.makeBasis(mesh, unitCell, group);
  
      TEST_ASSERT(basis.isInitialized());
      TEST_ASSERT(basis.isValid());
      TEST_ASSERT(eq(basis.nWave(), 4096));

      if (verbose() > 1) {
         std::ofstream out;
         openOutputFile("out/I_41_3_2", out);
         out << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(out);   
         basis.outputStars(out);   
         out.close();
      }

   }

};

TEST_BEGIN(BasisTest)
TEST_ADD(BasisTest, testConstructor)
TEST_ADD(BasisTest, testMake1DBasis_I)
TEST_ADD(BasisTest, testMake1DBasis_Ib)
TEST_ADD(BasisTest, testOutputWavesStars_1D)
TEST_ADD(BasisTest, testMake2DBasis_square_33)
TEST_ADD(BasisTest, testMake2DBasis_square_44)
TEST_ADD(BasisTest, testMake2DBasis_normTolerance)
TEST_ADD(BasisTest, testMake2DBasis_hex)
TEST_ADD(BasisTest, testMake3DBasis_I)
TEST_ADD(BasisTest, testMake3DBasis_I_m_3b_m)
TEST_ADD(BasisTest, testMake3DBasis_I_a_3b_d) 
TEST_ADD(BasisTest, testMake3DBasis_F_d_3b_m_2)
TEST_ADD(BasisTest, testMake3DBasis_F_d_3b_m_1)
TEST_ADD(BasisTest, testMake3DBasis_I_41_3_2)
TEST_END(BasisTest)

#endif
