#ifndef RP_SWEEP_FACTORY_H
#define RP_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>           // base class template
#include <rp/scft/sweep/Sweep.h>          // base class argument
#include <pscf/backend/TmplDeclare.h>    // declaration macros

#include <string>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Default Factory for subclasses of Sweep.
   *
   * \ingroup Rp_Scft_Sweep_Module
   */
   template <int D, class T>
   class SweepFactory : public Factory< Sweep<D,T> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      SweepFactory(System<D,T>& system);

      /**
      * Create an instance of a specified Sweep subclass.
      *
      * \param className  name of the Sweep subclass
      * \return  pointer to a new Sweep object
      */
      Sweep<D,T>* factory(std::string const & className) const;

      using Factory< Sweep<D,T> >::trySubfactories;
      using Factory< Sweep<D,T> >::readObjectOptional;

   private:

      // Pointer to parent system object.
      System<D,T>* systemPtr_;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(SweepFactory)

} // namespace Rp
} // namespace Pscf
#endif
