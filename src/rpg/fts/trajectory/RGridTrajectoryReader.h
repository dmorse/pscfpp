#ifndef RPG_RGRID_TRAJECTORY_READER_H
#define RPG_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/RGridTrajectoryReader.h> // direct base class 
#include <rpg/system/Types.h>                        // base class argument
#include <prdc/cuda/RField.h>                        // base class member
#include <rpg/fts/trajectory/TrajectoryReader.h>     // indirect base class 

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Trajectory file reader.
   *
   * \ingroup Rpg_Fts_Trajectory_Module
   */
   template <int D>
   class RGridTrajectoryReader 
    : public Rp::RGridTrajectoryReader<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      RGridTrajectoryReader<D>(Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RGridTrajectoryReader<1, Rpg::Types<1> >;
      extern template class RGridTrajectoryReader<2, Rpg::Types<2> >;
      extern template class RGridTrajectoryReader<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class RGridTrajectoryReader<1>;
      extern template class RGridTrajectoryReader<2>;
      extern template class RGridTrajectoryReader<3>;
   }
}
#endif
