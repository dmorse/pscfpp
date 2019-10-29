#ifndef FD1D_COMPOSITION_SWEEP_H
#define FD1D_COMPOSITION_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
   * Base class for solution along a line in parameter space.
   *
   * \ingroup Fd1d_Sweep_Module
   */
   class CompositionSweep : public Sweep
   {

   public:

      /**
      * Default constructor.
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
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Initialization at beginning sweep.
      */
      virtual void setup();

      /**
      * Iterate to solution.
      *
      * \param s path length coordinate, in [0,1]
      */
      virtual void setState(double s);

      /**
      * Output data to a running summary.
      *
      * \param out  output file, open for writing
      * \param i  integer index
      * \param s  value of path length parameter s
      */
      virtual void outputSummary(std::ostream& out, int i, double s);

   private:

      /**
      * Molecular volume fractions at beginning of sweep (s=0).
      */
      DArray<double> phi0_;

      /**
      * Change in molecule volume fractions over sweep s=[0,1].
      */
      DArray<double> dPhi_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
