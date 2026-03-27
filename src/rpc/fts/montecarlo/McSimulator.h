#ifndef RPC_MC_SIMULATOR_H
#define RPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McSimulator.h>      // base class template
#include <rpc/system/Types.h>                   // template argument
#include <rpc/fts/simulator/Simulator.h>        // indirect base class
#include <rpc/fts/montecarlo/McMoveManager.h>   // member of base class
#include <rpc/fts/analyzer/AnalyzerManager.h>   // member of base class

namespace Pscf {
namespace Rpc {

   /**
   * Monte Carlo simulator for PS-FTS.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::McSimulator , and
   * inherit their public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * \see \ref rp_McSimulator_page "Manual Page"
   * \see \ref rp_McSimulator_page (manual page)
   *
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class McSimulator : public Rp::McSimulator< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(System<D>& system);

      McSimulator() = delete;
      McSimulator(McSimulator<D> const & ) = delete;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McSimulator<1, Rpc::Types<1> >;
      extern template class McSimulator<2, Rpc::Types<2> >;
      extern template class McSimulator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class McSimulator<1>;
      extern template class McSimulator<2>;
      extern template class McSimulator<3>;
   }
}
#endif
