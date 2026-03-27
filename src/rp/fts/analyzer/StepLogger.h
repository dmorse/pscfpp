#ifndef RP_STEP_LOGGER_H
#define RP_STEP_LOGGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace Pscf {
namespace Rp {

   /**
   * Periodically write the step index to a log file.
   *
   * Specializations of this template are used as base classes for 
   * two closely analogous class templates, also named StepLogger, 
   * that are defined in the Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_StepLogger_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class StepLogger : public T::Analyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      StepLogger(typename T::Simulator& simulator, 
                 typename T::System& system);

      /**
      * Destructor.
      */
      virtual ~StepLogger()
      {}

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

   private:

      using AnalyzerT = typename T::Analyzer;

   };

}
}
#endif
