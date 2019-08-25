#ifndef PSSP_ITERATOR_TEST_H
#define PSSP_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/iterator/AmIterator.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/math/Constants.h>

#include <util/format/Dbl.h>

#include <pssp/System.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class IteratorTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor1()
   {
      AmIterator<1> AmItr;
      System<1> sys;
   }

   void testReadParam()
   {
      System<1> sys;

      #if 0
      char* argv[3];
      argv[0] = (char *) "myName";
      argv[1] = (char *) "-p";
      argv[2] = (char *) "in/param";
      optind = 1;
      sys.setOptions(3, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      TEST_ASSERT(eq(sys.iterator().epsilon(), 1e-3));
      TEST_ASSERT(sys.iterator().maxHist() == 3);
      TEST_ASSERT(sys.iterator().maxItr() == 40);
   }

   void testAllocate()
   {
      System<1> sys;


      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      //Allocate done automatically in system.tpp
   }

   void testComputeDeviation()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif
      
      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      // Set systemPtr_->wFields()
      MeshIterator<1> iter(sys.mesh().dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         sys.wFieldGrid(0)[iter.rank()] = cos(twoPi * 
                (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                 double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                 double(iter.position(2))/double(sys.mesh().dimension(2)))
              ); 
         sys.wFieldGrid(1)[iter.rank()] = sin(twoPi * 
                (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                 double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                 double(iter.position(2))/double(sys.mesh().dimension(2)))
              );
      }
      for (int i = 0; i < sys.mixture().nMonomer(); i++) {
         sys.fft().forwardTransform(sys.wFieldGrid(i), sys.wFieldDft(i));
         sys.basis().convertFieldDftToComponents(sys.wFieldDft(i),sys.wField(i));
      }

      //set systemPtr_->cField();
      DArray<double> wTemp;
      DArray<double> cTemp;
      double sum1;
      wTemp.allocate(sys.mixture().nMonomer());
      cTemp.allocate(sys.mixture().nMonomer());
      for ( int i = 0; i < sys.basis().nStar(); ++i) {
         sum1 = 0;
         for (int j = 0; j < sys.mixture().nMonomer(); ++j) {
            wTemp[j] = sys.wField(j)[i];
            sum1 += sys.wField(j)[i];
         }

         sum1 /= sys.mixture().nMonomer();


         sys.cField(0)[i] = sys.interaction().chiInverse(1,0) * 
                              (-sum1 + wTemp[1]);
         sys.cField(1)[i] = sys.interaction().chiInverse(1,0) *
                              (-sum1 + wTemp[0]);

      }
      
      //calculate deviation by hand();
      sys.iterator().computeDeviation();
      
      //print devHists_
      /*for(int i = 0; i < sys.mixture().nMonomer();i++){
         std::cout<<"THis is devfield of "<<i<<std::endl;
         for(int j = 0; j < sys.basis().nStar();j++){
            std::cout<<Dbl(sys.iterator().devHists_[0][i][j])<<std::endl;
         }
      }*/       
   }

   void testIsConverged1()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";

      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      // Set systemPtr_->wFields()
      MeshIterator<1> iter(sys.mesh().dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         sys.wFieldGrid(0)[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
         sys.wFieldGrid(1)[iter.rank()] = sin(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
      }
      for (int i = 0; i < sys.mixture().nMonomer(); i++) {
         sys.fft().forwardTransform(sys.wFieldGrid(i), sys.wFieldDft(i));
         sys.basis().convertFieldDftToComponents(sys.wFieldDft(i),sys.wField(i));
      }

      //set systemPtr_->cField();
      DArray<double> wTemp;
      DArray<double> cTemp;
      double sum1;
      wTemp.allocate(sys.mixture().nMonomer());
      cTemp.allocate(sys.mixture().nMonomer());
      for ( int i = 0; i < sys.basis().nStar(); ++i) {
         sum1 = 0;
         for (int j = 0; j < sys.mixture().nMonomer(); ++j) {
            wTemp[j] = sys.wField(j)[i];
            sum1 += sys.wField(j)[i];
         }

         sum1 /= sys.mixture().nMonomer();


         sys.cField(0)[i] = sys.interaction().chiInverse(1,0) * 
                              (-sum1 + wTemp[1]);
         sys.cField(1)[i] = sys.interaction().chiInverse(1,0) *
                              (-sum1 + wTemp[0]);

      }
      
      //calculate deviation
      sys.iterator().computeDeviation();
      
      //should be zero so converged automatically
      TEST_ASSERT(sys.iterator().isConverged());    
   }

   void testIsConverged2()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      // Set systemPtr_->wFields()
      MeshIterator<1> iter(sys.mesh().dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         sys.wFieldGrid(0)[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
         sys.wFieldGrid(1)[iter.rank()] = sin(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
      }
      for (int i = 0; i < sys.mixture().nMonomer(); i++) {
         sys.fft().forwardTransform(sys.wFieldGrid(i), sys.wFieldDft(i));
         sys.basis().convertFieldDftToComponents(sys.wFieldDft(i),sys.wField(i));
      }

      //set systemPtr_->cField();
      DArray<double> wTemp;
      DArray<double> cTemp;
      double xi;
      wTemp.allocate(sys.mixture().nMonomer());
      cTemp.allocate(sys.mixture().nMonomer());
      for ( int i = 0; i < sys.basis().nStar(); ++i) {
         for (int j = 0; j < sys.mixture().nMonomer(); ++j) {
            wTemp[j] = sys.wField(j)[i];
         }

         sys.interaction().computeC(wTemp, cTemp, xi);

         for (int j = 0; j < sys.mixture().nMonomer(); ++j) {
            sys.cField(j)[i] = cTemp[j];
         }
      }
      
      //calculate deviation by hand();
      sys.iterator().computeDeviation();
      
      //dev is not zero. check calculation by hand
      //print stuff
      sys.iterator().isConverged();
   }

   void testSolve()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      // Set systemPtr_->wFields()
      MeshIterator<1> iter(sys.mesh().dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         sys.wFieldGrid(0)[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
         sys.wFieldGrid(1)[iter.rank()] = sin(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
      }
      for (int i = 0; i < sys.mixture().nMonomer(); i++) {
         sys.fft().forwardTransform(sys.wFieldGrid(i), sys.wFieldDft(i));
         sys.basis().convertFieldDftToComponents(sys.wFieldDft(i),sys.wField(i));
      }

      int fail = sys.iterator().solve();
      std::cout<<std::endl;
      std::cout<<"Iterator fail? "<<fail<<std::endl;
      TEST_ASSERT(!fail);
   }

   void testOutput()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      MeshIterator<1> iter(sys.mesh().dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         sys.wFieldGrid(0)[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
         sys.wFieldGrid(1)[iter.rank()] = sin(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) );
      }

      for (int i = 0; i < sys.mixture().nMonomer(); i++) {
         sys.fft().forwardTransform(sys.wFieldGrid(i), sys.wFieldDft(i));
         sys.basis().convertFieldDftToComponents(sys.wFieldDft(i),sys.wField(i));
      }

      std::ofstream outFile;
      FileMaster().openOutputFile("omega.o", outFile);
      sys.writeFields(outFile, sys.wFields());
      outFile.close();
   }

   void testConvertFieldToGrid()
   {
      System<1> sys;

      #if 0
      char* argv[5];
      argv[0] = (char*) "diffName";
      argv[1] = (char*) "-p";
      argv[2] = (char*) "in/param";
      argv[3] = (char*) "-c";
      argv[4] = (char*) "in/command";
      optind = 1;
      sys.setOptions(5, argv);
      sys.readParam();
      #endif

      std::ifstream paramFile;
      openInputFile("in/param", paramFile);
      sys.readParam(paramFile);
      paramFile.close();

      std::ifstream inFile;
      FileMaster().openInputFile("omega.o", inFile);
      sys.readFields(inFile, sys.wFields());
      inFile.close();

      //convert to rgrid
      for (int i = 0; i < sys.mixture().nMonomer(); ++i) {
         sys.basis().convertFieldComponentsToDft(sys.wField(i), sys.wFieldDft(i));
         sys.fft().inverseTransform(sys.wFieldDft(i), sys.wFieldGrid(i));
      }

      std::ofstream outFile;
      FileMaster().openOutputFile("omegaGrid.o", outFile);
      sys.writeRFields(outFile, sys.wFieldGrids());
      outFile.close();

      //compare log with omegaGrid.o
      /*MeshIterator<3> iter(sys.mesh().dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         std::cout<<iter.rank()<<"    ";
         std::cout<<Dbl(cos(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) ));
         std::cout<<"    ";
         std::cout<< Dbl(sin(twoPi * 
                        (double(iter.position(0))/double(sys.mesh().dimension(0)) + 
                         double(iter.position(1))/double(sys.mesh().dimension(1)) + 
                         double(iter.position(2))/double(sys.mesh().dimension(2)) ) ));
         std::cout<<std::endl;
      }*/
   }

};


TEST_BEGIN(IteratorTest)
TEST_ADD(IteratorTest, testConstructor1)
TEST_ADD(IteratorTest, testReadParam)
TEST_ADD(IteratorTest, testAllocate)
TEST_ADD(IteratorTest, testComputeDeviation)
TEST_ADD(IteratorTest, testIsConverged1)
//TEST_ADD(IteratorTest, testIsConverged2)
TEST_ADD(IteratorTest, testSolve)
TEST_ADD(IteratorTest, testOutput)
TEST_ADD(IteratorTest, testConvertFieldToGrid)
TEST_END(IteratorTest)

#endif
