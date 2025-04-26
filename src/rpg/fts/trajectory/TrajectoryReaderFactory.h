#ifndef RPG_TRAJECTORY_READER_FACTORY_H
#define RPG_TRAJECTORY_READER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/trajectory/TrajectoryReader.h>
#include <string>


namespace Pscf {
namespace Rpg {
   
   template <int D> class System;
   
   using namespace Util;
   
   /**
   * Factory for subclasses of TrajectoryReader.
   *
   * \ingroup Rpg_Fts_Trajectory_Module
   */
   template <int D>
   class TrajectoryReaderFactory : public Factory< TrajectoryReader<D> > 
   {

   public:

      /// Constructor
      TrajectoryReaderFactory(System<D>& system);

      /**
      * Method to create any TrajectoryReader supplied with PSCF
      *
      * \param className name of the TrajectoryReader subclass
      * \return TrajectoryReader* pointer to new instance of className
      */
      TrajectoryReader<D>* factory(const std::string &className) const;
      
      using Factory< TrajectoryReader<D> >::trySubfactories;

   private:
      
      /// Pointer to the parent system.
      System<D>* sysPtr_;
      
   };

   #ifndef RPG_TRAJECTORY_READER_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class TrajectoryReaderFactory<1>;
   extern template class TrajectoryReaderFactory<2>;
   extern template class TrajectoryReaderFactory<3>;
   #endif

}
}
#endif
