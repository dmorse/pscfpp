/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <fd1d/System.h>
#include <fd1d/Domain.h>
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
    : ns_(0),
      homogeneousMode_(-1),
      baseFileName_(), 
      systemPtr_(0),
      mixturePtr_(0),
      domainPtr_(0),
      iteratorPtr_(0)
   {  setClassName("Sweep"); }

   Sweep::Sweep(System& system)
    : ns_(0),
      homogeneousMode_(-1),
      baseFileName_(), 
      systemPtr_(&system),
      mixturePtr_(&system.mixture()),
      domainPtr_(&system.domain()),
      iteratorPtr_(&system.iterator())
   {  setClassName("Sweep"); }

   Sweep::~Sweep()
   {}

   void Sweep::setSystem(System& system)
   {
      systemPtr_  = &system;
      mixturePtr_ = &(system.mixture());
      domainPtr_ = &(system.domain());
      iteratorPtr_ = &(system.iterator());
   }

   /*
   * Read parameters.
   */
   void Sweep::readParameters(std::istream& in)
   {
      read<int>(in, "ns", ns_);
      read<std::string>(in, "baseFileName", baseFileName_);
      readOptional<int>(in, "homogeneousMode", homogeneousMode_);
   }

   void Sweep::solve()
   {

      // Set Sweep object
      setup();

      // Open summary file
      std::ofstream outFile;
      std::string fileName = baseFileName_;
      fileName += ".out";
      system().fileMaster().openOutputFile(fileName, outFile);

      // Solve for initial state of sweep
      double s = 0.0;
      int i = 0;
      int error;
      error = iterator().solve();
      if (error) {
         UTIL_THROW("Failure to converge initial state of sweep");
      } else { 
         if (homogeneousMode_ >= 0) {
            system().computeHomogeneous(homogeneousMode_);
         }
         fileName = baseFileName_;
         fileName += ".";
         fileName += toString(i);
         outputSolution(fileName, s);
         outputSummary(outFile, i, s);
      }

      // Loop over states on path
      double ds = 1.0/double(ns_);
      double ds0 = ds;
      std::cout << "ns = " << ns_ << std::endl;
      std::cout << "ds = " << ds  << std::endl;
      bool finished = false;
      while (!finished) {
         error = 1;
         while (error) {

            std::cout << std::endl;
            std::cout << "Begin s = " << s << std::endl;
            setState(s+ds);
            error = iterator().solve();
            if (error) {
               ds *= 0.50;
               if (ds < 0.1*ds0) {
                  UTIL_THROW("Step size too small in sweep");
               }
            } else {
               if (homogeneousMode_ >= 0) {
                  system().computeHomogeneous(homogeneousMode_);
               }
               s += ds;
               ++i;
               fileName = baseFileName_;
               fileName += ".";
               fileName += toString(i);
               outputSolution(fileName, s);
               outputSummary(outFile, i, s);
            }
         }
         if (s + ds > 1.0001) {
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
      system().fileMaster().openOutputFile(outFileName, out);
      system().writeParam(out);
      out << std::endl;
      system().outputThermo(out);
      if (homogeneousMode_ >= 0) {
         system().outputHomogeneous(homogeneousMode_, out);
      }
      out.close();

      // Write concentration fields
      outFileName = fileName;
      outFileName += ".c";
      system().fileMaster().openOutputFile(outFileName, out);
      system().writeFields(out, system().cFields());
      out.close();
      
      // Write chemical potential fields
      outFileName = fileName;
      outFileName += ".w";
      system().fileMaster().openOutputFile(outFileName, out);
      system().writeFields(out, system().wFields());
      out.close();

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

} // namespace Fd1d
} // namespace Pscf
