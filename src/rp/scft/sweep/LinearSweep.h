#ifndef RP_LINEAR_SWEEP_H
#define RP_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"                   // base class
#include "SweepParameter.h"          // member
#include <util/containers/DArray.h>  // member
#include <util/global.h>
#include <iostream>

namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class BasisFieldState;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Sweep in which parameters vary linearly with sweep variable s.
   *
   * \see \ref scft_sweep_linear_sec "Manual page"
   * \ingroup Rp_Scft_Sweep_Module
   */
   template <int D, class T>
   class LinearSweep : public T::Sweep
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      LinearSweep(System<D,T>& system);

      ~LinearSweep() = default;

      /**
      * Read parameters from param file.
      *
      * \param in  parameter file input stream
      */
      void readParameters(std::istream& in);

      /**
      * Setup operation at the beginning of a sweep.
      *
      * Gets and stores initial values of individual swept parameters.
      */
      void setup();

      /**
      * Set state parameters before solving an SCFT problem.
      *
      * Called by SweepTmpl::sweep() for each state in a sweep.
      *
      * \param s  path length coordinate, in [0,1]
      */
      void setParameters(double s);

      /**
      * Output data to a running summary.
      *
      * \param out  output file, open for writing
      */
      void outputSummary(std::ostream& out);

   protected:

      // Inherited protected members
      using SweepT = typename T::Sweep;
      using SweepT::system;

   private:

      // Typename aliases
      using BasisFieldStateT = BasisFieldState<D,T>;
      using SweepTmplT = SweepTmpl< BasisFieldStateT >;

      /// Number of parameters being swept.
      int nParameter_;

      /// Array of SweepParameter objects.
      DArray< SweepParameter<D,T> > parameters_;

   };

}
}
#endif
