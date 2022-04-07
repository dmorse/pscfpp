/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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


   /*
   * Constructor.
   */
   CompositionSweep::CompositionSweep(System& system)
    : Sweep(system)
   {  setClassName("CompositionSweep"); }

   /*
   * Destructor.
   */
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
      Sweep::setup();
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         phi0_[i] = mixture().polymer(i).phi();
      }
   }

   /*
   * Set state for specified value of s.
   */
   void CompositionSweep::setParameters(double s)
   {
      double phi;
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         phi = phi0_[i] + s*dPhi_[i];
         mixture().polymer(i).setPhi(phi);
      }
   }

   void CompositionSweep::outputSummary(std::ostream& out)
   {
      double sNew = s(0);
      int i = nAccept() - 1;
      int np = mixture().nPolymer();
      if (homogeneousMode_ == -1) {
         out << Int(i) << Dbl(sNew) 
             << Dbl(system().fHelmholtz(), 20, 10)
             << Dbl(system().pressure(), 20, 10);
         for (int j = 0; j < np; ++j) {
            out << Dbl(mixture().polymer(j).phi(), 20, 10);
         }
         out << std::endl;
      } else {
         out << Int(i) << Dbl(sNew, 10);
         if (homogeneousMode_ == 0) {
            double dF = system().fHelmholtz() 
                      - system().homogeneous().fHelmholtz();
            out << Dbl(dF, 20, 10);
            for (int j = 0; j < np; ++j) {
               out << Dbl(mixture().polymer(j).phi(), 16);
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
            for (int j = 0; j < np; ++j) {
               out << Dbl(mixture().polymer(j).mu(), 16);
            }
            for (int j = 0; j < np; ++j) {
               out << Dbl(system().homogeneous().phi(j), 16);
            }
            double dV;
            for (int j = 0; j < np; ++j) {
               dV = mixture().polymer(j).phi()
                    - system().homogeneous().phi(j);
               dV *= V;
               out << Dbl(dV, 16);
            }
         }
         out << std::endl;
      }
   }

} // namespace Fd1d
} // namespace Pscf
