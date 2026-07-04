#ifndef RP_STEP_LOGGER_TPP
#define RP_STEP_LOGGER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/format/Int.h>
#include <util/misc/Log.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   StepLogger<D,T>::StepLogger(Simulator<D,T>& simulator, 
                             System<D,T>& system)
    : Analyzer<D,T>(simulator, system)
   {  ParamComposite::setClassName("StepLogger"); }

   /*
   * Read interval.
   */
   template <int D, class T>
   void StepLogger<D,T>::readParameters(std::istream& in)
   {  Analyzer<D,T>::readInterval(in); }

   /*
   * Periodically write the step index to a log file.
   */
   template <int D, class T>
   void StepLogger<D,T>::sample(long iStep)
   {
      if (Analyzer<D,T>::isAtInterval(iStep)) {
         Log::file() << "iStep  " << Int(iStep,10) << std::endl;
      }
   }

}
}
#endif
