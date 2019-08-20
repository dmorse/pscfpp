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
      //printEndl();

      Basis<3> sampleBasis;
      TEST_ASSERT(eq(sampleBasis.nWave(),0));
      TEST_ASSERT(eq(sampleBasis.nStar(),0));
      TEST_ASSERT(eq(sampleBasis.nBasis(),0));
   }



   void testMake1DBasis_1() {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<1> unitCell;
      std::ifstream in;
      openInputFile("in/Lamellar_8", in);
      in >> unitCell;

      // Make mesh object
      IntVec<1> d;
      in >> d;
      in.close();
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

   void testMake1DBasis_2() {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<1> unitCell;
      std::ifstream in;
      openInputFile("in/Lamellar_8", in);
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

      if (verbose() > 1) {
         basis.outputWaves(std::cout);
         basis.outputStars(std::cout);
      }

   }

   void testMake2DBasis_1()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/Square2D_33", in);
      in >> unitCell;

      // Make mesh object
      IntVec<2> d;
      in >> d;
      in.close();
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

   void testMake2DBasis_2()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/Square2D_44", in);
      in >> unitCell;

      // Make mesh object
      IntVec<2> d;
      in >> d;
      in.close();
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

      if (verbose() > 1) {
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
      }

   }

   void testMake2DBasis_3()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Read UnitCell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/hexagonal", in);
      in >> unitCell;
      in.close();

      // Make Mesh object
      IntVec<2> d;
      d[0] = 4;
      d[1] = 4;
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

      if (verbose() > 1) {
         std::cout << "nBasis = " << basis.nBasis() << std::endl;
         basis.outputWaves(std::cout);   
         basis.outputStars(std::cout);   
      }

   }

   void testFieldConversion2D()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<2> unitCell;
      std::ifstream in;
      openInputFile("in/Square2D_33", in);
      in >> unitCell;

      // Make mesh object
      IntVec<2> d;
      in >> d;
      in.close();
      Mesh<2> mesh(d);

      // Construct basis object
      Basis<2> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(basis.nWave(), 9));
      TEST_ASSERT(eq(basis.nStar(),9));
      TEST_ASSERT(basis.isValid());

      // Associate a grid wise field value
      RField<2> rField;
      RFieldDft<2> kField;
      RField<2> rFieldAlt;
      RFieldDft<2> kFieldAlt;
      DArray<double> components;
      FFT<2> fftMachine;

      components.allocate(mesh.size());
      rField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());

      fftMachine.setup(rField, kField);

      MeshIterator<2> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      double p0, p1, p11;
      for (iter.begin(); !iter.atEnd(); ++iter){  
         p0 = double(iter.position(0))/double(mesh.dimension(0));
         p0 *= twoPi;
         p1 = double(iter.position(1))/double(mesh.dimension(1));
         p1 *= twoPi;
         p11 = p0 + p1;
         rField[iter.rank()] = 1.3 
                             + 1.2*cos(p0) + 0.8*sin(p0) 
                             + 1.1*cos(p1) + 0.7*sin(p1)
                             + 1.4*cos(p11);
      }

      // std::cout<< "Forward Fourier transform" << std::endl;
      fftMachine.forwardTransform(rField, kField);

      #if 0
      std::cout << "Outputing kField values" << std::endl;
      for (int i = 0; i < kField.capacity(); i++) {
         std::cout  << kField[i][0]  << ' '
                    << kField[i][1]<<std::endl;
      }
      #endif
  
      // std::cout << "Converting DFT to basis functions" << std::endl;
      basis.convertFieldDftToComponents(kField, components);

      #if 0
      std::cout << "Outputing components" << std::endl;
      for( int i = 0; i < basis.nStar(); i++) {
         std::cout << Dbl(components[i]) 
                   << Dbl(components[i]/sqrt(2.0)) << std::endl;
      }
      #endif

      // std::cout << "Converting component fields to DFT" 
      //           << std::endl;
      basis.convertFieldComponentsToDft(components, kFieldAlt);

      #if 0
      std::cout << "Outputing kField values" << std::endl;
      for (int i = 0; i < kField.capacity(); i++) {
         std::cout  << kFieldAlt[i][0]  << ' '
                    << kFieldAlt[i][1]<<std::endl;
      }
      #endif
 
      for (int i = 0; i < kField.capacity(); i++) {
         TEST_ASSERT(eq(kField[i][0], kFieldAlt[i][0]));
         TEST_ASSERT(eq(kField[i][1], kFieldAlt[i][1]));
      }
 
      #if 0 
      fftMachine.inverseTransform(kFieldAlt, rFieldAlt);
      std::cout<<"Outputing rFieldAlt values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rFieldAlt[i]<<std::endl;
      }
      for ( int i = 0; i < mesh.size(); i++) {
         TEST_ASSERT(eq(rField[i], rFieldAlt[i]));
      }
      #endif

   }

   void testMake3DBasis()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell3D", in);
      in >> unitCell;

      // Make mesh object
      IntVec<3> d;
      in >> d;
      in.close();
      Mesh<3> mesh(d);

      // Construct basis object
      Basis<3> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(basis.nWave(), 27));
      TEST_ASSERT(eq(basis.nStar(),27));
      TEST_ASSERT(basis.isValid());
   }

   void testFieldConversion3D()
   {
      printMethod(TEST_FUNC);
      //printEndl();

      // Make unitcell
      UnitCell<3> unitCell;
      std::ifstream in;
      openInputFile("in/UnitCell3D", in);
      in >> unitCell;

      // Make mesh object
      IntVec<3> d;
      in >> d;
      in.close();
      Mesh<3> mesh(d);

      // Construct basis object
      Basis<3> basis;
      std::string spaceGroup = "I";
      basis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(basis.nWave(), 27));
      TEST_ASSERT(eq(basis.nStar(),27));

      // Associate a grid wise field value
      RField<3> rField;
      RFieldDft<3> kField;
      RField<3> rFieldAlt;
      RFieldDft<3> kFieldAlt;
      DArray<double> components;
      FFT<3> fftMachine;

      components.allocate(mesh.size());
      rField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());

      fftMachine.setup(rField, kField);

      MeshIterator<3> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      double p0, p1, p2, p11, p01, p02;
      for (iter.begin(); !iter.atEnd(); ++iter){  
         p0 = double(iter.position(0))/double(mesh.dimension(0));
         p0 *= twoPi;
         p1 = double(iter.position(1))/double(mesh.dimension(1));
         p1 *= twoPi;
         p2 = double(iter.position(2))/double(mesh.dimension(2));
         p2 *= twoPi;
         p01 = p0 + p1;
         p11 = p1 + p1;
         p02 = p0 + p2;
         p11 = p1 + p1;
         rField[iter.rank()] = 1.3 
                             + 1.2*cos(p0) + 0.8*sin(p0) 
                             + 1.1*cos(p1) + 0.7*sin(p1)
                             + 0.2*cos(p2) + 0.9*sin(p2)
                             + 1.4*cos(p11)
                             + 0.6*sin(p02);
      }

      // std::cout<< "Forward Fourier transform" << std::endl;
      fftMachine.forwardTransform(rField, kField);

      #if 0
      std::cout << "Outputing kField values" << std::endl;
      for (int i = 0; i < kField.capacity(); i++) {
         std::cout  << kField[i][0]  << ' '
                    << kField[i][1]<<std::endl;
      }
      #endif
  
      // std::cout << "Converting DFT to basis functions" << std::endl;
      basis.convertFieldDftToComponents(kField, components);

      #if 0
      std::cout << "Outputing components" << std::endl;
      for( int i = 0; i < basis.nStar(); i++) {
         std::cout << Dbl(components[i]) 
                   << Dbl(components[i]/sqrt(2.0)) << std::endl;
      }
      #endif

      // std::cout << "Converting component fields to DFT" 
      //           << std::endl;
      basis.convertFieldComponentsToDft(components, kFieldAlt);

      #if 0
      std::cout << "Outputing kField values" << std::endl;
      for (int i = 0; i < kField.capacity(); i++) {
         std::cout  << kFieldAlt[i][0]  << ' '
                    << kFieldAlt[i][1]<<std::endl;
      }
      #endif
 
      for (int i = 0; i < kField.capacity(); i++) {
         TEST_ASSERT(eq(kField[i][0], kFieldAlt[i][0]));
         TEST_ASSERT(eq(kField[i][1], kFieldAlt[i][1]));
      }
 
      #if 0 
      fftMachine.inverseTransform(kFieldAlt, rFieldAlt);
      std::cout<<"Outputing rFieldAlt values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rFieldAlt[i]<<std::endl;
      }

      for ( int i = 0; i < mesh.size(); i++) {
         TEST_ASSERT(eq(rField[i], rFieldAlt[i]));
      }
      #endif

   }

};

TEST_BEGIN(BasisTest)
TEST_ADD(BasisTest, testConstructor)
TEST_ADD(BasisTest, testMake1DBasis_1)
TEST_ADD(BasisTest, testMake1DBasis_2)
TEST_ADD(BasisTest, testMake2DBasis_1)
TEST_ADD(BasisTest, testMake2DBasis_2)
TEST_ADD(BasisTest, testMake2DBasis_3)
TEST_ADD(BasisTest, testFieldConversion2D)
TEST_ADD(BasisTest, testMake3DBasis)
TEST_ADD(BasisTest, testFieldConversion3D)
TEST_END(BasisTest)

#endif
