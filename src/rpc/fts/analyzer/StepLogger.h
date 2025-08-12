#ifndef RPC_STEP_LOGGER_H
#define RPC_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Periodically write snapshots to a trajectory file.
   *
   * \see rpc_StepLogger_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class StepLogger : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      StepLogger();

      /**
      * Destructor.
      */
      virtual ~StepLogger()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);

      using ParamComposite::setClassName;
      using Analyzer<D>::isAtInterval;

   };

   #ifndef RPC_STEP_LOGGER_TPP
   // Suppress implicit instantiation
   extern template class StepLogger<1>;
   extern template class StepLogger<2>;
   extern template class StepLogger<3>;
   #endif

}
}
#endif
