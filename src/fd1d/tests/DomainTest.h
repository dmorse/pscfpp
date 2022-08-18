#ifndef FD1D_DOMAIN_TEST_H
#define FD1D_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/domain/Domain.h>
#include <fd1d/domain/GeometryMode.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Fd1d;

class DomainTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Domain domain;
   }

   void testPlanarVolume()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMin = 0.7;
      double xMax = 2.5;
      int nx = 101;
      Domain domain;
      domain.setPlanarParameters(xMin, xMax, nx);
      TEST_ASSERT(eq(domain.volume(), xMax - xMin));
   }

   void testCylinderVolume()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMax = 2.0;
      int nx = 101;
      Domain domain;
      domain.setCylinderParameters(xMax, nx);
      double pi = Constants::Pi;
      double volume = pi*xMax*xMax;
      TEST_ASSERT(eq(domain.volume(), volume));

   }

   void testCylinderShellVolume()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMin = 1.5;
      double xMax = 2.0;
      int nx = 101;
      Domain domain;
      domain.setShellParameters(GeometryMode::Cylindrical, xMin, xMax, nx);
      double pi = Constants::Pi;
      double volume = pi*(xMax*xMax - xMin*xMin);
      TEST_ASSERT(eq(domain.volume(), volume));
   }

   void testSphereVolume()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMax = 2.0;
      int nx = 101;
      Domain domain;
      domain.setSphereParameters(xMax, nx);
      double pi = Constants::Pi;
      double volume = 4.0*pi*xMax*xMax*xMax/3.0;
      TEST_ASSERT(eq(domain.volume(), volume));

   }

   void testSphericalShellVolume()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMin = 1.5;
      double xMax = 2.0;
      int nx = 101;
      Domain domain;
      domain.setShellParameters(GeometryMode::Spherical, xMin, xMax, nx);
      double pi = Constants::Pi;
      double volume = 4.0*pi*(xMax*xMax*xMax - xMin*xMin*xMin)/3.0;
      TEST_ASSERT(eq(domain.volume(), volume));
   }

   void testCylindricalAverageUniform()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      int nx = 101;
      double xMax = 1.7;
      //double dx = xMax/double(nx-1);

      Domain domain;
      domain.setCylinderParameters(xMax, nx);

      DArray<double> f;
      f.allocate(nx);
      //double x;
      double A = 1.3;
      for (int i=0; i < nx; ++i) {
         //x = dx*double(i);
         f[i] = A;
      }
      if (verbose() > 0) {
         std::cout << "\n A       = " << A ;
         std::cout << "\n Average = " << domain.spatialAverage(f);
      }
      TEST_ASSERT(eq(domain.spatialAverage(f), A));
   }

   void testSphericalAverageUniform()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      int nx = 101;
      double xMax = 1.7;
      //double dx = xMax/double(nx-1);

      Domain domain;
      domain.setSphereParameters(xMax, nx);

      DArray<double> f;
      f.allocate(nx);
      //double x;
      double A = 1.3;
      for (int i=0; i < nx; ++i) {
         //x = dx*double(i);
         f[i] = A;
      }
      if (verbose() > 0) {
         std::cout << "\n A       = " << A ;
         std::cout << "\n Average = " << domain.spatialAverage(f);
      }
      TEST_ASSERT(eq(domain.spatialAverage(f), A));
   }

   void testCylindricalAverageLinear()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      int nx = 801;
      double xMax = 1.7;
      double dx = xMax/double(nx-1);

      Domain domain;
      domain.setCylinderParameters(xMax, nx);

      DArray<double> f;
      f.allocate(nx);
      double x;
      double B = 0.7;
      for (int i=0; i < nx; ++i) {
         x = dx*double(i);
         f[i] = B*x;
      }
      double computed  = domain.spatialAverage(f);
      double predicted = 2.0*B*xMax/3.0;
      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n computed   = " << computed;
         std::cout << "\n predicted  = " << predicted;
      }
      TEST_ASSERT(std::abs(computed - predicted) < 1.0E-4);
   }


   void testSphericalAverageLinear()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      int nx = 801;
      double xMax = 1.7;
      double dx = xMax/double(nx-1);

      Domain domain;
      domain.setSphereParameters(xMax, nx);

      DArray<double> f;
      f.allocate(nx);
      double x;
      double B = 0.7;
      for (int i=0; i < nx; ++i) {
         x = dx*double(i);
         f[i] = B*x;
      }
      double computed  = domain.spatialAverage(f);
      double predicted = 0.75*B*xMax;
      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n computed   = " << computed;
         std::cout << "\n predicted  = " << predicted;
      }
      TEST_ASSERT(std::abs(computed - predicted) < 1.0E-4);
   }

};

TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testConstructor)
TEST_ADD(DomainTest, testPlanarVolume)
TEST_ADD(DomainTest, testCylinderVolume)
TEST_ADD(DomainTest, testCylinderShellVolume)
TEST_ADD(DomainTest, testSphereVolume)
TEST_ADD(DomainTest, testSphericalShellVolume)
TEST_ADD(DomainTest, testCylindricalAverageUniform)
TEST_ADD(DomainTest, testSphericalAverageUniform)
TEST_ADD(DomainTest, testCylindricalAverageLinear)
TEST_ADD(DomainTest, testSphericalAverageLinear)
TEST_END(DomainTest)

#endif
