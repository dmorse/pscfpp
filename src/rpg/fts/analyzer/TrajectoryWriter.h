#ifndef RPG_TRAJECTORY_WRITER_H
#define RPG_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/TrajectoryWriter.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

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
   * \ingroup Rpg_Fts_Analyzer_Module
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
      TrajectoryWriter(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryWriter<1, Rpg::Types<1> >;
      extern template class TrajectoryWriter<2, Rpg::Types<2> >;
      extern template class TrajectoryWriter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class TrajectoryWriter<1>;
      extern template class TrajectoryWriter<2>;
      extern template class TrajectoryWriter<3>;
   }
}
#endif
