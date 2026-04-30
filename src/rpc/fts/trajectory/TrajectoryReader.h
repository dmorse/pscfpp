#ifndef RPC_TRAJECTORY_READER_H
#define RPC_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/TrajectoryReader.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   /**
   * Trajectory file reader (base class).
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::TrajectoryReader, and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \ingroup Rpc_Fts_Trajectory_Module
   */
   template <int D>
   class TrajectoryReader : public Rp::TrajectoryReader<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      TrajectoryReader<D>(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader<D>() = default;
   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryReader<1, Rpc::Types<1> >;
      extern template class TrajectoryReader<2, Rpc::Types<2> >;
      extern template class TrajectoryReader<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class TrajectoryReader<1>;
      extern template class TrajectoryReader<2>;
      extern template class TrajectoryReader<3>;
   }
}
#endif
