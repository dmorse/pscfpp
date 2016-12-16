/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CompositionSweep.h"
#include "System.h"
#include "Mixture.h"
#include "Domain.h"
#include "Iterator.h"

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
   * Read composition increment.
   */
   void CompositionSweep::readParameters(std::istream& in)
   {
      read<int>(in, "ns", ns_);
      int np = mixture().nPolymer();
      readDArray<double>(in, "dPhi", dPhi_, np);
   }

   /*
   * Setup operation at beginning sweep.
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
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         mixture().polymer(i).setPhi(phi0_[i] + s*dPhi_[i]);
      }
   }

   void CompositionSweep::solve()
   {

      // Solve for initial state of sweep
      int error;
      error = iterator().solve();
      if (error) {
         UTIL_THROW("Failure to converge initial state of sweep");
      }

      // Set Sweep object
      setup();

      // Loop over states on path
      double ds = 1.0/double(ns_);
      double ds0 = ds;
      double s = 0.0;
      bool finished = true;
      while (!finished) {
         error = 1;
         while (error && !finished) {
            setState(s+ds);
            error = iterator().solve();
            if (error) {
               ds *= 0.50;
               if (ds < 0.1*ds0) {
                  UTIL_THROW("Too small step size in sweeep");
               }
            } else {
               s += ds;
            }
         }
         if (s >= 1.0) {
            finished = true;
         }
      }
   }

} // namespace Fd1d
} // namespace Pscf
