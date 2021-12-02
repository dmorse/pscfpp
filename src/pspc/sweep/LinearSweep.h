#ifndef PSPC_LINEAR_SWEEP_H
#define PSPC_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"      // base class
#include "SweepParameter.h" // parameter class
#include <pspc/solvers/Block.h>
#include <pspc/solvers/Mixture.h>
#include <pspc/solvers/Polymer.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/global.h>
#include <iostream>
#include <algorithm>
#include <iomanip>

namespace Pscf {
namespace Pspc {

   // Forward declare classes and operators
   template <int D>
   class System;

   template <int D>
   class LinearSweep;

   template <int D>
   class SweepParameter;

   using namespace Util;

   /**
   * Base class for a sweep in parameter space where parameters change
   * linearly with the sweep variable. 
   * 
   * \ingroup Pspc_Sweep_Module
   */
   template <int D>
   class LinearSweep : public Sweep<D>
   {
   public:

      /** 
      * Constructor.
      * \param system parent System object
      */
      LinearSweep(System<D>& system);

      /**
      * Read parameters from param file.
      * 
      * \param in Input stream from param file.
      */
      void readParameters(std::istream& in);

      /**
      * Setup operation at the beginning of a sweep. Gets initial 
      * values of individual parameters.
      */
       void setup();

      /**
      * Set the state before an iteration. Called with each new iteration 
      * in SweepTempl::sweep()
      *
      * \param s path length coordinate, in [0,1]
      */    
      void setParameters(double s);

      /**
      * Output data to a running summary.
      *
      * \param out  output file, open for writing
      */
      void outputSummary(std::ostream& out);

   protected:

      using Sweep<D>::system;
      using Sweep<D>::hasSystem;
   
   private:

      /**
      * Number of parameters being swept. 
      */
      int nParameter_; 

      /**
      * Array of SweepParameter objects.
      */
      DArray< SweepParameter<D> > parameters_;

   // friends:
       
   };

   #ifndef PSPC_LINEAR_SWEEP_TPP
   // Suppress implicit instantiation
   extern template class LinearSweep<1>;
   extern template class LinearSweep<2>;
   extern template class LinearSweep<3>;
   #endif

}
}
#endif
