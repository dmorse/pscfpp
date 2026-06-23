#ifndef RPC_TRAJECTORY_WRITER_H
#define RPC_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/TrajectoryWriter.h> // direct base class template
#include <rpc/system/Types.h>                 // base class argument
#include "Analyzer.h"                         // indirect base class

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * Specializations of this template are derived from specializations of
   * the base class template Rp::TrajectoryWriter, and inherit their
   * entire public interface and almost all of their source code from
   * this base class.
   *
   * \see Rp::TrajectoryWriter
   * \see \ref rp_TrajectoryWriter_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class TrajectoryWriter : public Rp::TrajectoryWriter< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      TrajectoryWriter(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryWriter<1, Rpc::Types<1> >;
      extern template class TrajectoryWriter<2, Rpc::Types<2> >;
      extern template class TrajectoryWriter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class TrajectoryWriter<1>;
      extern template class TrajectoryWriter<2>;
      extern template class TrajectoryWriter<3>;
   }
}
#endif
