#ifndef RPC_RGRID_TRAJECTORY_READER_H
#define RPC_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/RGridTrajectoryReader.h> // direct base class 
#include <rpc/system/Types.h>                        // base class argument
#include <prdc/cpu/RField.h>                         // base class member
#include <rpc/fts/trajectory/TrajectoryReader.h>     // indirect base class 

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Trajectory file reader.
   *
   * \ingroup Rpc_Fts_Trajectory_Module
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
      RGridTrajectoryReader<D>(Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RGridTrajectoryReader<1, Rpc::Types<1> >;
      extern template class RGridTrajectoryReader<2, Rpc::Types<2> >;
      extern template class RGridTrajectoryReader<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class RGridTrajectoryReader<1>;
      extern template class RGridTrajectoryReader<2>;
      extern template class RGridTrajectoryReader<3>;
   }
}
#endif
