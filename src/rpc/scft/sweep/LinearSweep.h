#ifndef RPC_LINEAR_SWEEP_H
#define RPC_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"            // base class
#include "SweepParameter.h"   // member
#include <util/global.h>
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   using namespace Util;

   /**
   * Sweep in which parameters vary linearly with sweep variable s.
   *
   * See also: \ref scft_param_sweep_linear_sec "Parameter file format"
   *
   * \ingroup Rpc_Scft_Sweep_Module
   */
   template <int D>
   class LinearSweep : public Sweep<D>
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LinearSweep(System<D>& system);

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
      * Called by SweepTempl::sweep() for each new state in sweep.
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

      // Inherited public members
      using ParamComposite::read;
      using ParamComposite::readDArray;

   protected:

      // Inherited protected members
      using Sweep<D>::system;
      using Sweep<D>::hasSystem;
      using SweepTmpl< BasisFieldState<D> >::parameterTypes_;

   private:

      /// Number of parameters being swept.
      int nParameter_;

      /// Array of SweepParameter objects.
      DArray< SweepParameter<D> > parameters_;

   };

   #ifndef RPC_LINEAR_SWEEP_TPP
   // Suppress implicit instantiation
   extern template class LinearSweep<1>;
   extern template class LinearSweep<2>;
   extern template class LinearSweep<3>;
   #endif

}
}
#endif
