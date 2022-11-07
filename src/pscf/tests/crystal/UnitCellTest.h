#ifndef PSCF_UNIT_CELL_TEST_H
#define PSCF_UNIT_CELL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <util/math/Constants.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class UnitCellTest : public UnitTest 
{

public:

   void setUp()
   { setVerbose(0); }

   void tearDown()
   {}
 
   template <int D>
   bool isValidReciprocal(UnitCell<D> const & cell, bool verbose = false)
   {
      double sum;
      double twoPi = 2.0*Constants::Pi;
      bool isValid = true;
      int i, j, k;
      if (verbose) {
         Log::file() << std::endl;
      }
      for (i = 0; i < D; ++i ) {
         for (j = 0; j < D; ++j ) {
            sum = 0.0;
            for (k=0; k < D; ++k ) {
               sum += cell.rBasis(i)[k]*cell.kBasis(j)[k];  
            }
            sum = sum/twoPi;
            if (verbose) {
               Log::file() << Dbl(sum, 15, 5);
            }
            if (i == j) {
               sum -= 1.0;
            }
            if (std::abs(sum) > 1.0E-8) {
               isValid = false;
            }
         }
         if (verbose) {
            Log::file() << std::endl;
         }
      }
      return isValid;
   }

   template <int D>
   bool isValidDerivative(UnitCell<D> const & cell)
   {
      double sum;
      double nParams = cell.nParameter();
      int i, j, k, m;
      for (k = 0; k < nParams; ++k) {
         for (i = 0; i < D; ++i) {
            for (j = 0; j < D; ++j) {
               sum = 0.0;
               for (m = 0; m < D; ++m) {
                  sum += cell.drBasis(k, i, m)*cell.kBasis(j)[m];
                  sum += cell.dkBasis(k, j, m)*cell.rBasis(i)[m];
               }
            }
            if (std::abs(sum) > 1.0E-8) {
               return false;
            }
         }
      }
      return true;
   }

   void test1DLamellar() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<1> v;
      std::ifstream in;
      openInputFile("in/Lamellar", in);
      in >> v;
      double param = v.parameter(0);
      double twoPi = 2.0*Constants::Pi;

      TEST_ASSERT(eq(v.rBasis(0)[0], param));
      TEST_ASSERT(eq(v.kBasis(0)[0], twoPi/param));
      TEST_ASSERT(eq(v.volume(), param));
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test ksq function
      IntVec<1> x;
      int m = -3;
      x[0] = m;
      double xSq = v.ksq(x);
      double y = (twoPi*m)/param;
      TEST_ASSERT(eq(xSq, y*y));

      // Test assignment
      UnitCell<1> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      TEST_ASSERT(eq(u.rBasis(0)[0], param));
      TEST_ASSERT(eq(u.drBasis(0,0,0), 1.0));
      TEST_ASSERT(eq(u.kBasis(0)[0], twoPi/param));
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

      #if 1
      if (verbose() > 1) {
         std::cout.width(20);
         std::cout.precision(6);
         std::cout << v << std::endl ;
         std::cout << v.rBasis(0) << std::endl;
         std::cout << v.kBasis(0) << std::endl;
      }
      #endif

   }

   void test2DSquare() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Square", in);

      in >> v;
      TEST_ASSERT(v.nParameter() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double param = v.parameter(0);
      TEST_ASSERT(eq(v.volume(), param*param));

      #if 1
      if (verbose() > 0) {
         std::cout.width(20);
         std::cout.precision(6);
         std::cout << v << std::endl ;
   
         std::cout << "a(0) = " << v.rBasis(0) << std::endl;
         std::cout << "a(1) = " << v.rBasis(1) << std::endl;
         std::cout << "b(0) = " << v.kBasis(0) << std::endl;
         std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      }
      #endif

      double twoPi = 2.0*Constants::Pi;
      double b, dbb;
      int i, j, k;
      for (k = 0; k < v.nParameter(); ++k) {
         for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
               if (i == j) {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
               if (i == j) {
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }

      // Test assignment
      UnitCell<2> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 1);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test2DHexagonal() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Hexagonal2D", in);

      in >> v;
      TEST_ASSERT(v.nParameter() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      #if 1
      if (verbose() > 0) {
         std::cout.width(20);
         std::cout.precision(6);
         std::cout << v << std::endl ;
   
         std::cout << "a(0) = " << v.rBasis(0) << std::endl;
         std::cout << "a(1) = " << v.rBasis(1) << std::endl;
         std::cout << "b(0) = " << v.kBasis(0) << std::endl;
         std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      }
      #endif

      // Test volume
      double param = v.parameter(0);
      TEST_ASSERT(eq(v.volume(), param*param*0.5*sqrt(3.0)));

      IntVec<2> d;
      d[0] = 8;
      d[1] = 8;
      IntVec<2> x;
      x[0] = -4;
      x[1] = +5;
      IntVec<2> y;
     
      std::cout << "Before shift " << x << std::endl;
      y = shiftToMinimum(x, d, v);
      std::cout << "After shift  " << y << std::endl;
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == -3);
      y = shiftToMinimum(y, d, v);
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == -3);
      //std::cout << "After again  " << y << std::endl;

      // Test assignment
      UnitCell<2> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 1);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j;
      for (i = 0; i < 2; ++i) {
         for (j = 0; j < 2; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            TEST_ASSERT(eq(u.drBasis(0, i, j), v.drBasis(0, i, j)));
            TEST_ASSERT(eq(u.dkBasis(0, i, j), v.dkBasis(0, i, j)));
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test2DRectangular() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Rectangular", in);

      in >> v;
      TEST_ASSERT(v.nParameter() == 2);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double b = v.parameter(1);
      TEST_ASSERT(eq(v.volume(), a*b));

      // Test assignment
      UnitCell<2> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 2);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test2DRhombic() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Rhombic", in);

      in >> v;
      TEST_ASSERT(v.nParameter() == 2);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double alpha = v.parameter(1);
      TEST_ASSERT(eq(v.volume(), a*a*sin(alpha)));

      // Test assignment
      UnitCell<2> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 2);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j;
      for (i = 0; i < 2; ++i) {
         for (j = 0; j < 2; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            TEST_ASSERT(eq(u.drBasis(0, i, j), v.drBasis(0, i, j)));
            TEST_ASSERT(eq(u.dkBasis(0, i, j), v.dkBasis(0, i, j)));
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test2DOblique() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Oblique", in);

      in >> v;
      TEST_ASSERT(v.nParameter() == 3);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double b = v.parameter(1);
      double alpha = v.parameter(2);
      TEST_ASSERT(eq(v.volume(), a*b*sin(alpha)));

      // Test assignment
      UnitCell<2> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 3);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j;
      for (i = 0; i < 2; ++i) {
         for (j = 0; j < 2; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            TEST_ASSERT(eq(u.drBasis(0, i, j), v.drBasis(0, i, j)));
            TEST_ASSERT(eq(u.dkBasis(0, i, j), v.dkBasis(0, i, j)));
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test3DCubic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Cubic", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      TEST_ASSERT(eq(v.volume(), a*a*a));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "a(2) = " << v.rBasis(2) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      std::cout << "b(2) = " << v.kBasis(2) << std::endl;
      #endif

      #if 0
      double param, b, dbb;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (k = 0; k < v.nParameter(); ++k) {
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.parameter(i);
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.parameter(i);
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }
      #endif

      IntVec<3> d;
      d[0] = 8;
      d[1] = 8;
      d[2] = 8;
      IntVec<3> x;
      x[0] = -4;
      x[1] = +4;
      x[2] =  7;
      IntVec<3> y;
     
      //std::cout << "Before shift " << x << std::endl;
      y = shiftToMinimum(x, d, v);
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == 4);
      TEST_ASSERT(y[2] == -1);
      //std::cout << "After shift  " << y << std::endl;
      //y = shiftToMinimum(y, d, v);
      //std::cout << "After again  " << y << std::endl;

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 1);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test3DTetragonal() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Tetragonal", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 2);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double c = v.parameter(1);
      //double alpha = v.parameter(2);
      TEST_ASSERT(eq(v.volume(), a*a*c));

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 2);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j, k;
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test3DOrthorhombic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Orthorhombic", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 3);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double b = v.parameter(1);
      double c = v.parameter(2);
      TEST_ASSERT(eq(v.volume(), a*b*c));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "a(2) = " << v.rBasis(2) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      std::cout << "b(2) = " << v.kBasis(2) << std::endl;
      #endif

      double param, dbb;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (k = 0; k < v.nParameter(); ++k) {
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.parameter(i);
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.parameter(i);
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 3);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }


   void test3DHexagonal() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Hexagonal3D", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 2);

      bool isReciprocal = isValidReciprocal(v);          // Quiet test
      //bool isReciprocal = isValidReciprocal(v, true);  // Verbose test
      TEST_ASSERT(isReciprocal);
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double c = v.parameter(1);
      TEST_ASSERT(eq(v.volume(), a*a*c*0.5*sqrt(3.0)));

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 2);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j, k;
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }


   void test3DRhombohedral() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Rhombohedral", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 2);

      // Check rBasis, and interpretation of parameter a and beta
      double a = v.parameter(0);
      double beta = v.parameter(1);
      double sum00 = 0.0;
      double sum11 = 0.0;
      double sum22 = 0.0;
      double sum01 = 0.0;
      double sum12 = 0.0;
      for (int i=0; i < 3; ++i) {
          sum00 += v.rBasis(0)[i]*v.rBasis(0)[i];
          sum11 += v.rBasis(1)[i]*v.rBasis(1)[i];
          sum22 += v.rBasis(2)[i]*v.rBasis(2)[i];
          sum01 += v.rBasis(0)[i]*v.rBasis(1)[i];
          sum12 += v.rBasis(1)[i]*v.rBasis(2)[i];
      }
      TEST_ASSERT(eq(sum00, a*a));
      TEST_ASSERT(eq(sum11, a*a));
      TEST_ASSERT(eq(sum22, a*a));
      TEST_ASSERT(eq(sum01, a*a*cos(beta)));
      TEST_ASSERT(eq(sum12, a*a*cos(beta)));

      bool isReciprocal = isValidReciprocal(v);          // Quiet test
      //bool isReciprocal = isValidReciprocal(v, true);  // Verbose test
      TEST_ASSERT(isReciprocal);

      TEST_ASSERT(isValidDerivative(v));

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 2);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j, k;
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test3DMonoclinic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Monoclinic", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 4);

      bool isReciprocal = isValidReciprocal(v);          // Quiet test
      //bool isReciprocal = isValidReciprocal(v, true);  // Verbose test
      TEST_ASSERT(isReciprocal);
      TEST_ASSERT(isValidDerivative(v));

      // Test volume
      double a = v.parameter(0);
      double b = v.parameter(1);
      double c = v.parameter(2);
      double beta = v.parameter(3);
      TEST_ASSERT(eq(v.volume(), a*b*c*sin(beta)));

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 4);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j, k;
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

   void test3DTriclinic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Triclinic", in);
      in >> v;

      TEST_ASSERT(v.nParameter() == 6);

      bool isReciprocal = isValidReciprocal(v);          // Quiet test
      //bool isReciprocal = isValidReciprocal(v, true);  // Verbose test
      TEST_ASSERT(isReciprocal);

      TEST_ASSERT(isValidDerivative(v));

      // Test assignment
      UnitCell<3> u;
      u = v;
      TEST_ASSERT(u.lattice() == v.lattice());
      TEST_ASSERT(u.nParameter() == v.nParameter());
      TEST_ASSERT(u.nParameter() == 6);
      for (int i = 0; i < u.nParameter(); ++i) {
         TEST_ASSERT(eq(u.parameter(i), v.parameter(i)));
      }
      int i, j, k;
      for (i = 0; i < 3; ++i) {
         for (j = 0; j < 3; ++j) {
            TEST_ASSERT(eq(u.rBasis(i)[j], v.rBasis(i)[j]));
            TEST_ASSERT(eq(u.kBasis(i)[j], v.kBasis(i)[j]));
            for (k = 0; k < 3; ++k) {
               TEST_ASSERT(eq(u.drBasis(k,i,j), v.drBasis(k, i, j)));
               TEST_ASSERT(eq(u.dkBasis(k,i,j), v.dkBasis(k, i, j)));
            }
         }
      }
      TEST_ASSERT(isValidReciprocal(u));
      TEST_ASSERT(isValidDerivative(u));

   }

};

TEST_BEGIN(UnitCellTest)
TEST_ADD(UnitCellTest, test1DLamellar)
TEST_ADD(UnitCellTest, test2DSquare)
TEST_ADD(UnitCellTest, test2DHexagonal)
TEST_ADD(UnitCellTest, test2DRectangular)
TEST_ADD(UnitCellTest, test2DRhombic)
TEST_ADD(UnitCellTest, test2DOblique)
TEST_ADD(UnitCellTest, test3DCubic)
TEST_ADD(UnitCellTest, test3DTetragonal)
TEST_ADD(UnitCellTest, test3DOrthorhombic)
TEST_ADD(UnitCellTest, test3DHexagonal)
TEST_ADD(UnitCellTest, test3DRhombohedral)
TEST_ADD(UnitCellTest, test3DMonoclinic)
TEST_ADD(UnitCellTest, test3DTriclinic)
TEST_END(UnitCellTest)

#endif
