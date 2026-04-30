#ifndef RPC_SIMULATOR_H
#define RPC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>    // base class template
#include <rpc/system/Types.h>              // template argument
#include <rpc/fts/simulator/SimState.h>    // member
#include <prdc/cpu/RField.h>               // member (template arg)

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Field theoretic simulator (base class).
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::Simlator, and inherit
   * their entire public interface and almost all of their source code
   * from this base class.  See the documentation of this base class 
   * template for details. 
   *
   * Specializations of this template serve as as base classes for 
   * Rpc::BdSimulation<D> and Rpc::McSimulation<D>, for D=1, 2, and 3. 
   * For information about parameter file formats for subclasses, see:
   *
   * \see \ref rp_BdSimulator_page "MdSimulator manual page"
   * \see \ref rp_McSimulator_page "McSimulator manual page"
   * \ingroup Rpc_Fts_Simulator_Module
   */
   template <int D>
   class Simulator : public Pscf::Rp::Simulator<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      Simulator(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~Simulator() = default;

      // Prohibit copying and assignment.
      Simulator(Simulator<D> const &) = delete;
      Simulator<D>& operator = (Simulator<D> const &) = delete;

   private:

      /// Alias for direct base class.
      using RpSimulator = Pscf::Rp::Simulator<D, Types<D> >;

   };


} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Simulator<1, Rpc::Types<1> >;
      extern template class Simulator<2, Rpc::Types<2> >;
      extern template class Simulator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Simulator<1>;
      extern template class Simulator<2>;
      extern template class Simulator<3>;
   }
}
#endif
