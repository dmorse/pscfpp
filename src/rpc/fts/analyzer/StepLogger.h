#ifndef RPC_STEP_LOGGER_H
#define RPC_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                    // indirect base 
#include <rp/fts/analyzer/StepLogger.h>  // base class template
#include <rpc/system/Types.h>            // base class argument

namespace Pscf {
namespace Rpc {

   // Forward declarations
   template <int D> class Simulator;
   template <int D> class System;

   using namespace Util;

   /**
   * Periodically write the step index to a log file.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::StepLogger, and inherit their entire 
   * public interface and almost all of their source code from this 
   * base class. 
   *
   * \see Rp::StepLogger
   * \see \ref rp_StepLogger_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class StepLogger : public Rp::StepLogger< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      StepLogger(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class StepLogger< 1, Rpc::Types<1> >;
      extern template class StepLogger< 2, Rpc::Types<2> >;
      extern template class StepLogger< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class StepLogger<1>;
      extern template class StepLogger<2>;
      extern template class StepLogger<3>;
   }
}
#endif
