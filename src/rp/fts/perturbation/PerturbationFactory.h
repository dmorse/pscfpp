#ifndef RP_PERTURBATION_FACTORY_H
#define RP_PERTURBATION_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                 // base class template
#include <rp/fts/perturbation/Perturbation.h>   // base class argument

#include <pscf/backend/TmplDeclare.h>
#include <string>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Factory for subclasses of Perturbation.
   *
   * \ingroup Rp_Fts_Perturbation_Module
   */
   template <int D, class T>
   class PerturbationFactory
    : public Factory< Perturbation<D,T> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D,T> object
      */
      PerturbationFactory(Simulator<D,T>& simulator);

      /**
      * Method to create any Perturbation supplied with PSCF.
      *
      * \param className name of the Perturbation subclass
      * \return Perturbation* pointer to new instance of className
      */
      Perturbation<D,T>* factory(const std::string & className) const;

      using Factory< Perturbation<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Simulator<D,T>* simulatorPtr_;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(PerturbationFactory)

}
}
#endif
