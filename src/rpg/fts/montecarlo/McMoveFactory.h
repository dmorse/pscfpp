#ifndef RPG_MC_MOVE_FACTORY_H
#define RPG_MC_MOVE_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
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
   class McMoveFactory : public Factory< Rp::McMove<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Rp::McSimulator<D, Rpg::Types<D> > object
      */
      McMoveFactory(Rp::McSimulator<D, Rpg::Types<D> >& simulator);

      /**
      * Method to create any McMove supplied with PSCF.
      *
      * \param className name of the McMove subclass
      * \return McMove* pointer to new instance of className
      */
      Rp::McMove<D, Rpg::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::McMove<D, Rpg::Types<D> > >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Rp::McSimulator<D, Rpg::Types<D> >* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class McMoveFactory<1>;
   extern template class McMoveFactory<2>;
   extern template class McMoveFactory<3>;

}
}
#endif
