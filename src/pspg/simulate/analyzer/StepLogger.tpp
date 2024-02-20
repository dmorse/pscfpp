#ifndef PSPG_STEP_LOGGER_TPP
#define PSPG_STEP_LOGGER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StepLogger.h"
#include "Analyzer.h"
#include <util/format/Int.h>
#include <util/misc/Log.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   StepLogger<D>::StepLogger()
    : Analyzer<D>()
   {  setClassName("StepLogger"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void StepLogger<D>::readParameters(std::istream& in) 
   {  Analyzer<D>::readInterval(in); }

   /*
   * Periodically write a frame to file
   */
   template <int D>
   void StepLogger<D>::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Log::file() << "iStep  " << Int(iStep,10) << std::endl;
      }
   }

}
}
#endif
