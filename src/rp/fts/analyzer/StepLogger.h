#ifndef RP_STEP_LOGGER_H
#define RP_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h>

#include <pscf/backend/TmplDeclare.h>
#include <iostream>

namespace Pscf {
namespace Rp {

   /**
   * Periodically write the step index to a log file.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_StepLogger_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class StepLogger : public Analyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      StepLogger(Simulator<D,T>& simulator, 
                 System<D,T>& system);

      /**
      * Destructor.
      */
      ~StepLogger() = default;

      /**
      * Read interval.
      *
      * \param in input parameter file
      */
      void readParameters(std::istream& in) override;

      /**
      * Write the step index to a log file.
      *
      * \param iStep  step index
      */
      void sample(long iStep) override;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(StepLogger)

}
}
#endif

