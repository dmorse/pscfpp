/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MuSweep.h"
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

   MuSweep::MuSweep()
    : Sweep()
   {  setClassName("MuSweep"); }

   MuSweep::MuSweep(System& system)
    : Sweep(system)
   {  setClassName("MuSweep"); }

   MuSweep::~MuSweep()
   {}

   /*
   * Read parameters.
   */
   void MuSweep::readParameters(std::istream& in)
   {
      // Read ns, baseFileName and (optionally) homogeneousMode
      Sweep::readParameters(in);

      int np = mixture().nPolymer();
      dMu_.allocate(np);
      mu0_.allocate(np);
      readDArray<double>(in, "dMu", dMu_, np);
   }

   /*
   * Initialization at beginning sweep. Set mu0 to current composition.
   */
   void MuSweep::setup()
   {
      Sweep::setup();
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         UTIL_CHECK(mixture().polymer(i).ensemble() == Species::Open);
         mu0_[i] = mixture().polymer(i).mu();
      }
   }

   /*
   * Set state for specified value of s.
   */
   void MuSweep::setParameters(double s)
   {
      double mu;
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         mu = mu0_[i] + s*dMu_[i];
         mixture().polymer(i).setMu(mu);
      }
   }

   void MuSweep::outputSummary(std::ostream& out) 
   {
      int i = nAccept() - 1;
      int np = mixture().nPolymer();
      double sNew = s(0);
      if (homogeneousMode_ == -1) {
         out << Dbl(sNew) 
             << Dbl(system().fHelmholtz(), 20, 10)
             << Dbl(system().pressure(), 20, 10);
         for (i = 0; i < np - 1; ++i) {
            out << Dbl(mixture().polymer(i).phi(), 20, 10);
         }
         out << std::endl;
      } else {
         out << Dbl(sNew, 10);
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
   }

} // namespace Fd1d
} // namespace Pscf
