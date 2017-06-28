#ifndef PSSP_BASIS_TEST_H
#define PSSP_BASIS_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

#include <pssp/basis/Basis.h>
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
   {}

   void tearDown()
   {}

   void testConstructor()
   {
      Basis<3> sampleBasis;
      TEST_ASSERT(eq(sampleBasis.nWave(),0));
      TEST_ASSERT(eq(sampleBasis.nStar(),0));
      TEST_ASSERT(eq(sampleBasis.nBasis(),0));
   }

   void testMakeBasis2D(){
      Basis<2> sampleBasis;

      UnitCell<2> unitCell;
      //make unitcell
      std::ifstream in;
      openInputFile("in/UnitCell2D", in);
      in >> unitCell;
      IntVec<2> d;
      in >> d;
      in.close();

      //make mesh object
      Mesh<2> mesh(d);
      std::string spaceGroup = "I";
      sampleBasis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(sampleBasis.nWave(), 9));
      TEST_ASSERT(eq(sampleBasis.nStar(),9));
      TEST_ASSERT(eq(sampleBasis.nBasis(), 9));

      //shotgun test for unsorted + I spacegroup
      MeshIterator<2> itr(mesh.dimensions());
      itr.begin();
      for(int i = 0; i < mesh.size(); i++) {

         int size = sampleBasis.star(i).size;
         int bId = sampleBasis.star(i).beginId;
         int eId = sampleBasis.star(i).endId;
         TEST_ASSERT(eq(size,1));
         TEST_ASSERT(eq(bId,eId));

         TEST_ASSERT(!sampleBasis.star(i).cancel); 

         Basis<2>::Wave test1 = sampleBasis.wave(i);
         Basis<2>::Wave test2 = sampleBasis.wave(itr.position());
         ++itr;
         TEST_ASSERT(eq(test1.starId, test2.starId));
      }
      

      /*std::cout<<"starID\t"<<"Index\t"<<"sqNorm\t"<<"Implicit\t"
               <<"InvertFlag\t"<<"Coeff\t"<<std::endl;
      for(int i = 0; i < mesh.size(); i++) {
         std::cout<<sampleBasis.wave(i).starId<<"\t";
         std::cout<<sampleBasis.star(i).waveBz[0]<<"\t";
         std::cout<<sampleBasis.star(i).waveBz[1]<<"\t";
         std::cout<<sampleBasis.wave(i).sqNorm<<"\t";
         std::cout<<sampleBasis.wave(i).implicit<<"\t";
         std::cout<<sampleBasis.star(i).invertFlag<<"\t";
         std::cout<<sampleBasis.wave(i).coeff<<"\t";
         std::cout<<sampleBasis.wave(i).indicesDft[0]<<"\t";
         std::cout<<sampleBasis.wave(i).indicesDft[1]<<"\t";           
         std::cout<<std::endl;
      }*/

   }

   void testFieldConversion2D()
   {
      Basis<2> sampleBasis;

      UnitCell<2> unitCell;
      //make unitcell
      std::ifstream in;
      openInputFile("in/UnitCell2D", in);
      in >> unitCell;
      IntVec<2> d;
      in >> d;
      in.close();

      //make mesh object
      Mesh<2> mesh(d);
      std::string spaceGroup = "I";

      //construct basis
      sampleBasis.makeBasis(mesh, unitCell, spaceGroup);

      //associate a grid wise field value
      RField<2> rField;
      RFieldDft<2> kField;
      RField<2> rFieldAlt;
      RFieldDft<2> kFieldAlt;
      DArray<double> components;
      FFT<2> fftMachine;

      components.allocate(mesh.size());

      /*//There is a bug here that would break the code if allocated in this order
      rField.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());*/

      //This is for some reason fine
      rField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());

      fftMachine.setup(rField, kField);

      MeshIterator<2> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         //std::cout<<iter.position(0)<<iter.position(1)<<iter.position(2)<<std::endl;
         rField[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) ) ); 
                        //+ 
                         //double(iter.position(2))/double(mesh.dimension(2)) ) );
      }

      /*std::cout<<"Outputing rField values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rField[i]<<std::endl;
      }*/

      fftMachine.forwardTransform(rField, kField);

      /*std::cout<<"Outputing kField values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<kField[i][0]<<' '<<kField[i][1]<<std::endl;
      }*/
      
      sampleBasis.convertFieldDftToComponents(kField, components);

      /*std::cout<<"Converting DFT fields to basis functions"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<Dbl(components[i])<<std::endl;
      }*/

      sampleBasis.convertFieldComponentsToDFT(components, kFieldAlt);

      /*std::cout<<"Converting component fields to DFT"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<kFieldAlt[i][0]<<' '<<kFieldAlt[i][1]<<std::endl;
      }*/
      /*for ( int i = 0; i < mesh.size(); i++) {
         //std::cout<<kField[i][0]<<' '<<kField[i][1]<<std::endl;
         TEST_ASSERT(eq(kField[i][0],kFieldAlt[i][0]));
         TEST_ASSERT(eq(kField[i][1],kFieldAlt[i][1]));
      }*/


      fftMachine.inverseTransform(kFieldAlt, rFieldAlt);
      /*std::cout<<"Outputing rFieldAlt values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rFieldAlt[i]<<std::endl;
      }*/

      for ( int i = 0; i < mesh.size(); i++) {
         TEST_ASSERT(eq(rField[i], rFieldAlt[i]));
      }


   }

   void testMakeBasis(){
      Basis<3> sampleBasis;

      UnitCell<3> unitCell;
      //make unitcell
      std::ifstream in;
      openInputFile("in/UnitCell", in);
      in >> unitCell;
      IntVec<3> d;
      in >> d;
      in.close();

      //make mesh object
      Mesh<3> mesh(d);
      std::string spaceGroup = "I";
      sampleBasis.makeBasis(mesh, unitCell, spaceGroup);
      
      TEST_ASSERT(eq(sampleBasis.nWave(), 27));
      TEST_ASSERT(eq(sampleBasis.nStar(),27));
      TEST_ASSERT(eq(sampleBasis.nBasis(), 27));

      //shotgun test for unsorted + I spacegroup
      MeshIterator<3> itr(mesh.dimensions());
      itr.begin();
      for(int i = 0; i < mesh.size(); i++) {

         int size = sampleBasis.star(i).size;
         int bId = sampleBasis.star(i).beginId;
         int eId = sampleBasis.star(i).endId;
         TEST_ASSERT(eq(size,1));
         TEST_ASSERT(eq(bId,eId));

         TEST_ASSERT(!sampleBasis.star(i).cancel); 

         Basis<3>::Wave test1 = sampleBasis.wave(i);
         Basis<3>::Wave test2 = sampleBasis.wave(itr.position());
         ++itr;
         TEST_ASSERT(eq(test1.starId, test2.starId));
      }
      

      /*std::cout<<"starID\t"<<"Index\t"<<"sqNorm\t"<<"Implicit\t"
               <<"InvertFlag\t"<<"Coeff\t"<<std::endl;
      for(int i = 0; i < mesh.size(); i++) {
         std::cout<<sampleBasis.wave(i).starId<<"\t";
         std::cout<<sampleBasis.star(i).waveBz[0]<<"\t";
         std::cout<<sampleBasis.star(i).waveBz[1]<<"\t";
         std::cout<<sampleBasis.star(i).waveBz[2]<<"\t";
         std::cout<<sampleBasis.wave(i).sqNorm<<"\t";
         std::cout<<sampleBasis.wave(i).implicit<<"\t";
         std::cout<<sampleBasis.star(i).invertFlag<<"\t";
         std::cout<<sampleBasis.wave(i).coeff<<"\t";
         std::cout<<sampleBasis.wave(i).indicesDft[0]<<"\t";
         std::cout<<sampleBasis.wave(i).indicesDft[1]<<"\t";    
         std::cout<<sampleBasis.wave(i).indicesDft[2]<<"\t";          
         std::cout<<std::endl;
      }*/

   }

   void testFieldConversion()
   {
      Basis<3> sampleBasis;

      UnitCell<3> unitCell;
      //make unitcell
      std::ifstream in;
      openInputFile("in/UnitCell", in);
      in >> unitCell;
      IntVec<3> d;
      in >> d;
      in.close();

      //make mesh object
      Mesh<3> mesh(d);
      std::string spaceGroup = "I";

      //construct basis
      sampleBasis.makeBasis(mesh, unitCell, spaceGroup);

      //associate a grid wise field value
      RField<3> rField;
      RFieldDft<3> kField;
      RField<3> rFieldAlt;
      RFieldDft<3> kFieldAlt;
      DArray<double> components;
      FFT<3> fftMachine;

      components.allocate(mesh.size());

      /*//There is a bug here that would break the code if allocated in this order
      rField.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());*/

      //This is for some reason fine
      rField.allocate(mesh.dimensions());
      rFieldAlt.allocate(mesh.dimensions());
      kField.allocate(mesh.dimensions());
      kFieldAlt.allocate(mesh.dimensions());

      fftMachine.setup(rField, kField);

      MeshIterator<3> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         //std::cout<<iter.position(0)<<iter.position(1)<<iter.position(2)<<std::endl;
         rField[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) + 
                         double(iter.position(2))/double(mesh.dimension(2)) ) );
      }

     /* std::cout<<"Outputing rField values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rField[i]<<std::endl;
      }*/

      fftMachine.forwardTransform(rField, kField);

      /*std::cout<<"Outputing kField values"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<kField[i][0]<<' '<<kField[i][1]<<std::endl;
      }*/
      
      sampleBasis.convertFieldDftToComponents(kField, components);

      /*std::cout<<"Converting DFT fields to basis functions"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<Dbl(components[i])<<std::endl;
      }*/

      sampleBasis.convertFieldComponentsToDFT(components, kFieldAlt);

      /*std::cout<<"Converting component fields to DFT"<<std::endl;
      for( int i = 0; i < mesh.size(); i++) {
         std::cout<<kFieldAlt[i][0]<<' '<<kFieldAlt[i][1]<<std::endl;
      }*/
      /*for ( int i = 0; i < mesh.size(); i++) {
         //std::cout<<kField[i][0]<<' '<<kField[i][1]<<std::endl;
         TEST_ASSERT(eq(kField[i][0],kFieldAlt[i][0]));
         TEST_ASSERT(eq(kField[i][1],kFieldAlt[i][1]));
      }*/

      //std::cout<<"Outputing rFieldAlt values"<<std::endl;
      fftMachine.inverseTransform(kFieldAlt, rFieldAlt);
      /*for( int i = 0; i < mesh.size(); i++) {
         std::cout<<rFieldAlt[i]<<std::endl;
      }*/

      for ( int i = 0; i < mesh.size(); i++) {
         TEST_ASSERT(eq(rField[i], rFieldAlt[i]));
      }


   }



};

TEST_BEGIN(BasisTest)
TEST_ADD(BasisTest, testConstructor)
TEST_ADD(BasisTest, testMakeBasis2D)
TEST_ADD(BasisTest, testFieldConversion2D)
TEST_ADD(BasisTest, testMakeBasis)
TEST_ADD(BasisTest, testFieldConversion)
TEST_END(BasisTest)

#endif
