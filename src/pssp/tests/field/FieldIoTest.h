#ifndef PSSP_FIELD_IO_TEST_H
#define PSSP_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/RField.h>
#include <pssp/field/RFieldDft.h>
#include <pssp/field/FFT.h>
#include <pscf/crystal/Basis.h>

#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class FieldIoTest : public UnitTest 
{

public:

   void setUp()
   {
      setVerbose(2);
   }

   void tearDown()
   {}

   void testFieldConversion2D()
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

   void testFieldConversion3D()
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

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testFieldConversion2D)
TEST_ADD(FieldIoTest, testFieldConversion3D)
TEST_END(FieldIoTest)

#endif
