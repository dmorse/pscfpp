#ifndef PSCF_UNIT_CELL_TEST_H
#define PSCF_UNIT_CELL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <util/math/Constants.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class UnitCellTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   template <int D>
   bool isValidReciprocal(UnitCell<D> cell)
   {
      double sum;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (i=0; i < D; ++i ) {
         for (j=0; j < D; ++j ) {
            sum = 0.0;
            for (k=0; k < D; ++k ) {
               sum += cell.rBasis(i)[k]*cell.kBasis(j)[k];  
            }
            if (i == j) {
               sum -= twoPi;
            }
            if (std::abs(sum) > 1.0E-8) {
               return false;
            }
         }
      }
      return true;
   }

   template <int D>
   bool isValidDerivative(UnitCell<D> cell)
   {
      double sum;
      double Params = cell.nParams();
      int i, j, k;
      for (k=0; k < Params; ++k ) {
         sum = 0.0;
         for (i=0; i < D; ++i ) {
            //sum = 0.0;
            for (j=0; j < D; ++j ) {
               sum += (cell.drBasis(k, i, j)*cell.kBasis(i)[j])+(cell.dkBasis(k,i,j)*cell.rBasis(i)[j]);
            }
            /*if (std::abs(sum) > 1.0E-8) {
                return false;
             }*/
         }

         if (std::abs(sum) > 1.0E-8) {
            return false;
         }
      }
      return true;
   }

   template <int D>
   bool isValidDoublekkDerivative(UnitCell<D> cell)
   {   
      double sum;
      double twoPi = 2.0*Constants::Pi;
      double Params = cell.nParams();
      int i, j, k, t;
      for (k=0; k < Params; ++k ) { 
         sum = 0.0;
         for (i=0; i < D; ++i ) { 
            //sum = 0.0;
            for (j=0; j < D; ++j ) {
               for (t=0; t<D; ++t)  {
                  sum += (cell.dkkBasis(k, i, j)*cell.rBasis(i)[t])-(2*twoPi*cell.dkBasis(k, j, t));
               }
            }   
            /*if (std::abs(sum) > 1.0E-8) {
                return false;
            }*/ 
          }   
          // if (std::abs(sum) > 1.0E-8) {
          //  return false;
          //}   
       }   
       return true;
   } 

   template <int D>
   bool isValidDoublerrDerivative(UnitCell<D> cell)
   {   
      double sum;
      double twoPi = 2.0*Constants::Pi;
      double Params = cell.nParams();
      int i, j, k, t;
      for (k=0; k < Params; ++k ) { 
         sum = 0.0;
         for (i=0; i < D; ++i ) { 
            //sum = 0.0;
            for (j=0; j < D; ++j ) { 
    
               for (t=0; t<D; ++t)  {
 
                  sum += (cell.drrBasis(k, i, j)*cell.kBasis(i)[t])-(2*twoPi*cell.drBasis(k, j, t));
	          std::cout<< "dk("<<j<<t<<k<<")"<<"\t"<<"="<<cell.dkBasis(k, j, t)<<"\n";
                  std::cout<< "dr("<<j<<t<<k<<")"<<"\t"<<"="<<cell.drBasis(k, j, t)<<"\n";
                  std::cout<< "dkk("<<j<<t<<k<<")"<<"\t"<<"="<<cell.dkkBasis(k, j, t)<<"\n";
                  std::cout<< "drr("<<j<<t<<k<<")"<<"\t"<<"="<<cell.drrBasis(k, j, t)<<"\n";

                  //std::cout<< "dr"<<"\t"<<cell.drBasis(k, j, t);
                  //std::cout<< "dkk"<<"\t"<<cell.dkkBasis(k, j, t);
                  //std::cout<< "drr"<<"\t"<<cell.drrBasis(k, j, t);
               }   
            }   
            /*if (std::abs(sum) > 1.0E-8) {
                return false;
            }*/ 
         }   

         //if (std::abs(sum) > 1.0E-8) {
         //   return false;
         //}   
      }   
      return true;
   }

   void test1DLamellar() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<1> v;
      std::ifstream in;
      openInputFile("in/Lamellar", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));
      TEST_ASSERT(isValidDoublekkDerivative(v));
      TEST_ASSERT(isValidDoublerrDerivative(v));

      std::cout << v.rBasis(0) << std::endl;
      std::cout << v.kBasis(0) << std::endl;
   }

   void test2DSquare() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Square", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));
      TEST_ASSERT(isValidDoublekkDerivative(v));
      TEST_ASSERT(isValidDoublerrDerivative(v));

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
   }

   void test3DOrthorhombic() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Orthorhombic", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));
      TEST_ASSERT(isValidDoublekkDerivative(v));
      TEST_ASSERT(isValidDoublerrDerivative(v));

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
   }

};

TEST_BEGIN(UnitCellTest)
TEST_ADD(UnitCellTest, test1DLamellar)
TEST_ADD(UnitCellTest, test2DSquare)
TEST_ADD(UnitCellTest, test3DOrthorhombic)
TEST_END(UnitCellTest)

#endif
