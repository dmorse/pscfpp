#ifndef FD1D_LENGTH_SWEEP_H
#define FD1D_LENGTH_SWEEP_H

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
   class LengthSweep : public Sweep
   {

   public:

      /**
      * Default constructor.
      */
      LengthSweep();

      /**
      * Constructor.
      */
      LengthSweep(System& system);

      /**
      * Destructor.
      */
      ~LengthSweep();

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
      virtual void setParameters(double s);

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
      *  Id of polymer that is being swept
      */
      int polymerId_;

      /**
      *  Id of block in the polymer that is being swept
      */
      int blockId_;

      /**
      * Length homopolymer at start of sweep (s=0).
      */
      double length0_;

      /**
      * Change in kuhn length over sweep s=[0,1].
      */
      double dLength_;


   };

} // namespace Fd1d
} // namespace Pscf
#endif
