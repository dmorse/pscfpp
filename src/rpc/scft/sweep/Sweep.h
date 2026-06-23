#ifndef RPC_SWEEP_H
#define RPC_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/Sweep.h>             // base class template
#include <rpc/system/Types.h>                // base class argument
#include <rpc/scft/sweep/BasisFieldState.h>  // indirect base argument

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Solve a sequence of SCFT problems along a line in parameter space.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base class template Rp::Sweep, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::Sweep
   * \see \ref scft_sweep_page "Manual page"
   * \ingroup Rpc_Scft_Sweep_Module
   */
   template <int D>
   class Sweep : public Rp::Sweep<D, Types<D> >
   {

   public:

      /**
      * Default Constructor.
      */
      Sweep();

      /**
      * Constructor, creates assocation with parent system.
      *
      * \param system  parent system
      */
      Sweep(Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Sweep() = default;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   extern template class SweepTmpl< Rp::BasisFieldState<1, Rpc::Types<1> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<2, Rpc::Types<2> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<3, Rpc::Types<3> > >;
   namespace Rp {
      extern template class Sweep<1, Rpc::Types<1> >;
      extern template class Sweep<2, Rpc::Types<2> >;
      extern template class Sweep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Sweep<1>;
      extern template class Sweep<2>;
      extern template class Sweep<3>;
   }
}
#endif
