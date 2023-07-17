#ifndef PSPC_TRAJECTORY_READER_FACTORY_H
#define PSPC_TRAJECTORY_READER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <pspc/simulate/trajectory/TrajectoryReader.h>
#include <pspc/System.h>
#include <pspc/simulate/McSimulator.h>
#include <string>


namespace Pscf {
namespace Pspc {

   using namespace Util;
   
   /**
   * Factory for subclasses of TrajectoryReader.
   *
   * \ingroup Pspc_Analyzer_Module
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

}
}
#endif
