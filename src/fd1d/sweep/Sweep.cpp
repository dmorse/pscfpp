/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <fd1d/System.h>
#include <fd1d/domain/Domain.h>
#include <fd1d/solvers/Mixture.h>
#include <fd1d/iterator/Iterator.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   Sweep::Sweep()
    : SystemAccess(),
      ns_(0),
      homogeneousMode_(-1),
      baseFileName_(),
      comparison_(),
      fieldIo_()
   {  setClassName("Sweep"); }

   Sweep::Sweep(System& system)
    : SystemAccess(system),
      ns_(0),
      homogeneousMode_(-1),
      baseFileName_(),
      comparison_(system),
      fieldIo_(system)
   {  setClassName("Sweep"); }

   Sweep::~Sweep()
   {}

   void Sweep::setSystem(System& system)
   {
      SystemAccess::setSystem(system);
      comparison_.setSystem(system);
      fieldIo_.setSystem(system);
   }

   /*
   * Read parameters.
   */
   void Sweep::readParameters(std::istream& in)
   {
      read<int>(in, "ns", ns_);
      read<std::string>(in, "baseFileName", baseFileName_);
      homogeneousMode_ = -1; // default value
      readOptional<int>(in, "homogeneousMode", homogeneousMode_);
   }

   void Sweep::solve()
   {

      int nm = mixture().nMonomer();
      int nx = domain().nx();
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(nx > 0);

      // Allocate memory for solutions
      if (!wFields0_.isAllocated()) {
         wFields0_.allocate(nm);
         for (int i = 0; i < nm; ++i) {
            wFields0_[i].allocate(nx);
         }
      }
      if (!wFields1_.isAllocated()) {
         wFields1_.allocate(nm);
         for (int i = 0; i < nm; ++i) {
            wFields1_[i].allocate(nx);
         }
      }

      // Compute and output ds
      double ds = 1.0/double(ns_);
      double ds0 = ds;
      std::cout << std::endl;
      std::cout << "ns = " << ns_ << std::endl;
      std::cout << "ds = " << ds  << std::endl;

      // Set Sweep object
      setup();

      // Open summary file
      std::ofstream outFile;
      std::string fileName = baseFileName_;
      fileName += "log";
      fileMaster().openOutputFile(fileName, outFile);

      // Solve for initial state of sweep
      double s = 0.0;
      int i = 0;
      int error;
      std::cout << std::endl;
      std::cout << "Begin s = " << s << std::endl;
      bool isContinuation = false; // False on first step
      error = system().iterator().solve(isContinuation);
      if (error) {
         UTIL_THROW("Failure to converge initial state of sweep");
      } else {
         assignFields(wFields0_, wFields());
         if (homogeneousMode_ >= 0) {
            comparison_.compute(homogeneousMode_);
         }
         fileName = baseFileName_;
         fileName += toString(i);
         outputSolution(fileName, s);
         outputSummary(outFile, i, s);
      }

      // Loop over states on path
      bool finished = false;   // Are we finished with the loop?
      int  nPrev = 0;          // Number of previous solutions stored
      double s1 = 0.0;         // Value of s at previous solution (if any)
      while (!finished) {
         error = 1;
         while (error) {

            std::cout << std::endl;
            std::cout << "Attempt s = " << s + ds << std::endl;

            // Setup guess for fields
            if (nPrev == 0) {
               std::cout << "Zeroth order continuation" << std::endl;
            } else {
               // std::cout << "1st order continuation" << std::endl;
               double f1 = ds/(s - s1);
               double f0 = 1.0 + f1;
               // std::cout << " ds = " << ds << std::endl;
               // std::cout << "  s = " << s  << ", s1 = " << s1 << std::endl;
               // std::cout << " f1 = " << f1 << ", f0 = " << f0 << std::endl;
               int i, j;
               for (i = 0; i < nm; ++i) {
                  for (j = 0; j < nx; ++j) {
                     wFields()[i][j] = f0*wFields0_[i][j] - f1*wFields1_[i][j];
                  }
               }
            }

            // Attempt solution
            setState(s+ds);
            isContinuation = true;
            error = system().iterator().solve(isContinuation);

            if (error) {

               // Upon failure, reset to fields from last converged solution
               assignFields(wFields(), wFields0_);

               // Decrease ds by half
               ds *= 0.50;
               if (ds < 0.2*ds0) {
                  UTIL_THROW("Step size too small in sweep");
               }

            } else {

               // Upon success, save new field
               s1 = s;
               nPrev = 1;
               assignFields(wFields1_, wFields0_);
               assignFields(wFields0_, wFields());

               // Compare to homogeneous reference system
               if (homogeneousMode_ >= 0) {
                  comparison_.compute(homogeneousMode_);
               }

               // Update s and output
               s += ds;
               ++i;
               fileName = baseFileName_;
               fileName += toString(i);
               outputSolution(fileName, s);
               outputSummary(outFile, i, s);

            }
         }
         if (s + ds > 1.0000001) {
            finished = true;
         }
      }
   }

   void Sweep::outputSolution(std::string const & fileName, double s)
   {
      std::ofstream out;
      std::string outFileName;

      // Write parameter file, with thermodynamic properties at end
      outFileName = fileName;
      outFileName += ".prm";
      fileMaster().openOutputFile(outFileName, out);
      system().writeParam(out);
      out << std::endl;
      system().outputThermo(out);
      if (homogeneousMode_ >= 0) {
         comparison_.output(homogeneousMode_, out);
      }
      out.close();

      // Write concentration fields
      outFileName = fileName;
      outFileName += ".c";
      fieldIo_.writeFields(cFields(), outFileName);

      // Write chemical potential fields
      outFileName = fileName;
      outFileName += ".w";
      fieldIo_.writeFields(wFields(), outFileName);
   }

   void Sweep::outputSummary(std::ostream& out, int i, double s)
   {
      if (homogeneousMode_ == -1) {
      out << Int(i,5) << Dbl(s)
          << Dbl(system().fHelmholtz(),16)
          << Dbl(system().pressure(),16)
          << std::endl;
      } else {
         out << Int(i,5) << Dbl(s)
             << Dbl(system().fHelmholtz(),16)
             << Dbl(system().pressure(),16);
         if (homogeneousMode_ == 0) {
            double dF = system().fHelmholtz()
                      - system().homogeneous().fHelmholtz();
            out << Dbl(dF, 16);
         } else {
            double dP = system().pressure()
                      - system().homogeneous().pressure();
            double dOmega = -1.0*dP*domain().volume();
            out << Dbl(dOmega, 16);
         }
         out << std::endl;
      }
   }

   void Sweep::assignFields(DArray<System::Field>& lhs,
                            DArray<System::Field> const & rhs) const
   {

      int nm = mixture().nMonomer();
      int nx = domain().nx();

      UTIL_CHECK(lhs.capacity() == nm);
      UTIL_CHECK(rhs.capacity() == nm);
      int i, j;
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(rhs[i].isAllocated());
         UTIL_CHECK(rhs[i].capacity() == nx);
         UTIL_CHECK(lhs[i].isAllocated());
         UTIL_CHECK(lhs[i].capacity() == nx);
         for (j = 0; j < nx; ++j) {
            lhs[i][j] = rhs[i][j];
         }
      }
   }

} // namespace Fd1d
} // namespace Pscf
