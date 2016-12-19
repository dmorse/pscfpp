/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CompositionSweep.h"
#include <fd1d/System.h>
#include <fd1d/Domain.h>
#include <fd1d/solvers/Mixture.h>
#include <fd1d/iterator/Iterator.h>

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
      Sweep::readParameters(in);
      // read<int>(in, "ns", ns_);
      // read<std::string>(in, "baseFileName", baseFileName_);

      int np = mixture().nPolymer();
      dPhi_.allocate(np);
      phi0_.allocate(np);
      readDArray<double>(in, "dPhi", dPhi_, np);
   }

   /*
   * Initialization at beginning sweep.
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

} // namespace Fd1d
} // namespace Pscf
