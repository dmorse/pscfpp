#ifndef RP_RAMP_FACTORY_H
#define RP_RAMP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>   // base class template
#include <rp/fts/ramp/Ramp.h>     // base class argument

#include <pscf/backend/TmplDeclare.h>
#include <string>

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class Simulator;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Factory for subclasses of Ramp.
   *
   * \ingroup Rp_Fts_Ramp_Module
   */
   template <int D, class T>
   class RampFactory : public Factory< Ramp<D,T> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D, T > object
      */
      RampFactory(Simulator<D,T>& simulator);

      /**
      * Method to create any Ramp supplied with PSCF.
      *
      * \param className name of the Ramp subclass
      * \return Ramp* pointer to new instance of className
      */
      Ramp<D,T>* factory(const std::string & className) const;

      using Factory< Ramp<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Simulator<D,T>* simulatorPtr_;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(RampFactory)

} // namespace Rp
} // namespace Pscf
#endif
