#ifndef RP_TRAJECTORY_READER_FACTORY_H
#define RP_TRAJECTORY_READER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <rp/fts/trajectory/TrajectoryReader.h>  // template argument

#include <pscf/backends/TmplDeclare.h>
#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Factory for subclasses of TrajectoryReader.
   *
   * \ingroup Rp_Fts_Trajectory_Module
   */
   template <int D, class T>
   class TrajectoryReaderFactory 
     : public Factory< TrajectoryReader<D,T> >
   {

   public:

      /*
      * Constructor.
      *
      * \param system  parent system object
      */
      TrajectoryReaderFactory(System<D,T>& system);

      /**
      * Create and instance of a specified subclass of TrajectoryReader.
      *
      * \param className  name of the TrajectoryReader subclass
      * \return TrajectoryReader  pointer to new instance of className
      */
      TrajectoryReader<D,T>* factory(const std::string &className) const;

      using Factory< TrajectoryReader<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent system.
      System<D,T>* sysPtr_;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(TrajectoryReaderFactory)

}
}
#endif
