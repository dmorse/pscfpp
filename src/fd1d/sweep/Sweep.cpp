/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
//#include <pscf/sweep/SweepTmpl.h>
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
    : Base(),
      SystemAccess(),
      homogeneousMode_(-1),
      comparison_(),
      fieldIo_()
   {  setClassName("Sweep"); }

   Sweep::Sweep(System& system)
    : Base(),
      SystemAccess(system),
      homogeneousMode_(-1),
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
      //SweepTmpl<Sweep::State>::readParameters(in);
      Base::readParameters(in);
      homogeneousMode_ = -1; // default value
      readOptional<int>(in, "homogeneousMode", homogeneousMode_);
   }

   /**
   * Setup operation at beginning sweep.
   *
   * Must call initializeHistory.
   */
   void Sweep::setup() 
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

      initializeHistory(wFields0_, wFields1_);

      // Open log summary file
      std::string fileName = baseFileName_;
      fileName += "log";
      fileMaster().openOutputFile(fileName, logFile_);

   };

   /**
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   void Sweep::setParameters(double s) {};

   /**
   * Create guess for adjustable variables by continuation.
   */
   void Sweep::setGuess(double sNew) 
   {
      int nm = mixture().nMonomer();
      int nx = domain().nx();
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(nx > 0);
      if (nAccept() == 1) {
         std::cout << "Zeroth order continuation" << std::endl;
      } else {
         double ds = sNew - s(0);
         double f1 = ds/(sNew - s(1));
         double f0 = 1.0 + f1;
         int i, j;
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wFields()[i][j] = f0*state(0)[i][j] - f1*state(1)[i][j];
            }
         }
      }
   };

   /**
   * Call current iterator to solve SCFT problem.
   *
   * Return 0 for sucessful solution, 1 on failure to converge.
   */
   int Sweep::solve(bool isContinuation) 
   {  return system().iterator().solve(isContinuation); };

   /**
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   void Sweep::reset() 
   {  assignFields(wFields(), state(0)); };

   /**
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current 
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   void Sweep::getSolution() 
   {
      // Assign current wFields to state(0)
      assignFields(state(0), wFields());

      if (homogeneousMode_ >= 0) {
         comparison_.compute(homogeneousMode_);
      }

      // Output solution
      int i = nAccept() - 1;
      std::string fileName = baseFileName_;
      fileName += toString(i);
      outputSolution(fileName);

      // Output brief summary to log file
      outputSummary(logFile_);
   };

   void Sweep::outputSolution(std::string const & fileName)
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

   void Sweep::outputSummary(std::ostream& out)
   {
      int i = nAccept() - 1;
      double sNew = s(0);
      out << Int(i,5) << Dbl(sNew)
          << Dbl(system().fHelmholtz(),16)
          << Dbl(system().pressure(),16);
      #if 0
      if (homogeneousMode_ != -1) {
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
      }
      #endif
      out << std::endl;
   }

   void Sweep::cleanup()
   {  logFile_.close(); }

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
