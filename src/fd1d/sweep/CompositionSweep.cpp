/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CompositionSweep.h"
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

   CompositionSweep::CompositionSweep()
    : Sweep()
   {  setClassName("CompositionSweep"); }

   CompositionSweep::CompositionSweep(System& system)
    : Sweep(system)
   {  setClassName("CompositionSweep"); }

   CompositionSweep::~CompositionSweep()
   {}

   /*
   * Read parameters.
   */
   void CompositionSweep::readParameters(std::istream& in)
   {
      // Read ns, baseFileName and (optionally) homogeneousMode
      Sweep::readParameters(in);

      int np = mixture().nPolymer();
      dPhi_.allocate(np);
      phi0_.allocate(np);
      readDArray<double>(in, "dPhi", dPhi_, np);
   }

   /*
   * Initialization at beginning sweep. Set phi0 to current composition.
   */
   void CompositionSweep::setup()
   {
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         phi0_[i] = mixture().polymer(i).phi();
      }
   }

   /*
   * Set state for specified value of s.
   */
   void CompositionSweep::setState(double s)
   {
      double phi;
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         phi = phi0_[i] + s*dPhi_[i];
         mixture().polymer(i).setPhi(phi);
      }
   }

   void CompositionSweep::outputSummary(std::ostream& out, int i, double s) 
   {
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
            double V = domain().volume();
            double fExV = fEx*V;
            double pExV = pEx*V;
            out << Dbl(fExV, 20, 10);
            out << Dbl(pExV, 20, 10);
            for (int i = 0; i < np; ++i) {
               out << Dbl(system().mixture().polymer(i).mu(), 16);
            }
            double dV;
            for (int i = 0; i < np - 1; ++i) {
               dV = system().mixture().polymer(i).phi()
                  - system().homogeneous().phi(i);
               dV *= V;
               out << Dbl(dV, 16);
            }
         }
         out << std::endl;
      }
   }

} // namespace Fd1d
} // namespace Pscf
