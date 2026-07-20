#ifndef RP_MC_MOVE_FACTORY_H
#define RP_MC_MOVE_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rp/fts/montecarlo/McMove.h>

#include <string>

namespace Pscf {
namespace Rp {

   // Forward declaration
   template <int D, class T> class McSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of McMove.
   *
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class McMoveFactory : public Factory< McMove<D,T> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator<D,T> object
      */
      McMoveFactory(McSimulator<D,T>& simulator);

      /**
      * Method to create any McMove supplied with PSCF.
      *
      * \param className name of the McMove subclass
      * \return McMove* pointer to new instance of className
      */
      McMove<D,T>* factory(const std::string &className) const;

      using Factory< McMove<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      McSimulator<D,T>* simulatorPtr_;

   };

}
}
#endif
