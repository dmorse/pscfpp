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
      * Destructor.
      */
      ~CompositionSweep();

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve(){};

   };

} // namespace Fd1d
} // namespace Pscf
#endif
