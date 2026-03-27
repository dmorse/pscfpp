#ifndef RPC_BD_SIMULATOR_H
#define RPC_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdSimulator.h>       // base class template
#include <rpc/system/Types.h>                  // template argument
#include <rpc/fts/simulator/Simulator.h>       // indirect base class
#include <rpc/fts/analyzer/AnalyzerManager.h>  // member of base class

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   /**
   * Brownian dynamics simulator for PS-FTS.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::BdSimulator, and 
   * inherit their entire public interface and almost all of their source 
   * code from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see \ref rp_BdSimulator_page "Manual Page"
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class BdSimulator : public Rp::BdSimulator<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      BdSimulator(System<D>& system);

      BdSimulator() = delete;
      BdSimulator(BdSimulator<D> const &) = delete;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdSimulator<1, Rpc::Types<1> >;
      extern template class BdSimulator<2, Rpc::Types<2> >;
      extern template class BdSimulator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BdSimulator<1>;
      extern template class BdSimulator<2>;
      extern template class BdSimulator<3>;
   }
}
#endif
