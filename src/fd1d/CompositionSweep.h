#ifndef FD1D_COMPOSITION_SWEEP_H
#define FD1D_COMPOSITION_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"          // base class
#include <util/global.h>                  

namespace Pscf {
namespace Fd1d
{

   class System;
   using namespace Util;

   /**
   * Base class for classes solve along a line in parameter space.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class CompositionSweep : public Sweep
   {

   public:

      /**
      * Constructor.
      */
      CompositionSweep();

      /**
      * Constructor.
      */
      CompositionSweep(System& system);

      /**
      * Destructor.
      */
      ~CompositionSweep();

      /**
      * Read parameters.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Setup operation at beginning sweep.
      */
      virtual void setup();

      /**
      * Iterate to solution.
      *
      * \param s path length coordinate, in [0,1]
      */
      virtual void setState(double s);

      /**
      * Iterate to solution.
      */
      virtual void solve();

   private:

      DArray<double> phi0_;

      DArray<double> dPhi_;

      int ns_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
