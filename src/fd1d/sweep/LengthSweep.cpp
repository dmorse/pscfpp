/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LengthSweep.h"
#include <fd1d/System.h>
#include <fd1d/domain/Domain.h>
#include <fd1d/solvers/Mixture.h>
#include <fd1d/iterator/Iterator.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   LengthSweep::LengthSweep()
    : Sweep()
   {  setClassName("LengthSweep"); }

   LengthSweep::LengthSweep(System& system)
    : Sweep(system)
   {  setClassName("LengthSweep"); }

   LengthSweep::~LengthSweep()
   {}

   /*
   * Read parameters.
   */
   void LengthSweep::readParameters(std::istream& in)
   {
      // Read ns, baseFileName and (optionally) homogeneousMode
      Sweep::readParameters(in);

      read<int>(in,"polymerId",polymerId_);
      read<int>(in,"blockId",blockId_);
      read<double>(in, "dLength", dLength_);
   }

   /*
   * Initialization at beginning sweep. Set L0 to current block.
   */
   void LengthSweep::setup()
   {
      length0_ = mixture().polymer(polymerId_).block(blockId_).length();
   }

   /*
   * Set state for specified value of s.
   */
   void LengthSweep::setState(double s)
   {
      double L = length0_ + s*dLength_;
      mixture().polymer(polymerId_).block(blockId_).setLength(L);
      //std::cout <<" s        = " <<  s 
      //          << " Polymer = " << mixture().polymer(polymerId_).length() 
      //          <<" Block    = " << mixture().polymer(polymerId_).block(blockId_).length() 
      //          << '\n';
   }

   void LengthSweep::outputSummary(std::ostream& out, int i, double s) 
   {   
      #if 0
      int np = mixture().nPolymer();
      if (homogeneousMode_ == -1) {
         out << Dbl(s) 
             << Dbl(system().fHelmholtz(), 20, 10)
             << Dbl(system().pressure(), 20, 10);
         for (i = 0; i < np - 1; ++i) {
            out << Dbl(mixture().polymer(i).phi(), 20, 10);
         }
         out << std::endl;
      } else {
         out << Dbl(s,10);
         if (homogeneousMode_ == 0) {
            double dF = system().fHelmholtz() 
                      - system().homogeneous().fHelmholtz();
            out << Dbl(dF, 20, 10);
            for (int i = 0; i < np - 1; ++i) {
               out << Dbl(mixture().polymer(i).phi(), 16);
            }
         } else {
            double fEx = system().fHelmholtz() 
                       - system().homogeneous().fHelmholtz();
            double pEx = system().pressure() 
                       - system().homogeneous().pressure();
            double V = domain().volume()/mixture().vMonomer();
            double fExV = fEx*V;
            double pExV = pEx*V;
            out << Dbl(fExV, 20, 10);
            out << Dbl(pExV, 20, 10);
            for (int i = 0; i < np; ++i) {
               out << Dbl(mixture().polymer(i).mu(), 16);
            }
            for (int i = 0; i < np - 1; ++i) {
               out << Dbl(system().homogeneous().phi(i), 16);
            }
            double dV;
            for (int i = 0; i < np - 1; ++i) {
               dV = mixture().polymer(i).phi()
                  - system().homogeneous().phi(i);
               dV *= V;
               out << Dbl(dV, 16);
            }
         }
         out << std::endl;
      }
      #endif
   }

} // namespace Fd1d
} // namespace Pscf
