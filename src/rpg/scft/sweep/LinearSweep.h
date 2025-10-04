#ifndef RPG_LINEAR_SWEEP_H
#define RPG_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"            // base class
#include "SweepParameter.h"   // member 
#include <util/global.h>
#include <iostream>

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Base class for a sweep in parameter space where parameters change
   * linearly with the sweep variable. 
   * 
   * \sa \ref scft_sweep_linear_sec "Parameter file format"
   *
   * \ingroup Rpg_Scft_Sweep_Module
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
      using SweepTmpl< BasisFieldState<D> >::parameterTypes_;
   
   private:
      /// Number of parameters being swept. 
      int nParameter_; 

      /// Array of SweepParameter objects.
      DArray< SweepParameter<D> > parameters_;
       
   };

   // Explicit instantiation declarations
   extern template class LinearSweep<1>;
   extern template class LinearSweep<2>;
   extern template class LinearSweep<3>;

}
}
#endif
