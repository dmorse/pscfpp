#ifndef RPG_STEP_LOGGER_H
#define RPG_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/StepLogger.h>  // base class template
#include <rpg/system/Types.h>            // base class argument
#include <rpg/fts/analyzer/Analyzer.h>   // indirect base

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Periodically write the step index to a log file.
   *
   * Specializations of this template are derived from specializations of
   * the base class template Rp::StepLogger, and inherit their entire
   * public interface and almost all of their source code from this base
   * class.
   *
   * \see Rp::StepLogger
   * \see \ref rp_StepLogger_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
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
      StepLogger(Rp::Simulator<D, Rpg::Types<D> >& simulator, Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class StepLogger< 1, Rpg::Types<1> >;
      extern template class StepLogger< 2, Rpg::Types<2> >;
      extern template class StepLogger< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class StepLogger<1>;
      extern template class StepLogger<2>;
      extern template class StepLogger<3>;
   }
}
#endif
