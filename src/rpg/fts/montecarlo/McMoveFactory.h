#ifndef RPG_MC_MOVE_FACTORY_H
#define RPG_MC_MOVE_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/montecarlo/McMove.h>

#include <string>

namespace Pscf {
namespace Rpg {

   template <int D> class McSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of McMove.
   *
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class McMoveFactory : public Factory< McMove<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator<D> object
      */
      McMoveFactory(McSimulator<D>& simulator);

      /**
      * Method to create any McMove supplied with PSCF.
      *
      * \param className name of the McMove subclass
      * \return McMove* pointer to new instance of className
      */
      McMove<D>* factory(const std::string &className) const;

      using Factory< McMove<D> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      McSimulator<D>* simulatorPtr_;

   };

   #ifndef RPG_MC_MOVE_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class McMoveFactory<1>;
   extern template class McMoveFactory<2>;
   extern template class McMoveFactory<3>;
   #endif

}
}
#endif
