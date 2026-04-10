#ifndef RPG_MIXTURE_H
#define RPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Mixture.h>    // base class template
#include <rpg/system/Types.h>      // base class argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of base class template Rp::Mixture, and
   * inherit their public interface and almost all of their source 
   * code from this base class.
   *
   * \see Rp::Mixture
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Mixture : public Rp::Mixture< D, Types<D> >
   {

   public:

      /// Base class type aliases
      using RpMixtureT = typename Rp::Mixture<D, Types<D> >;
      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::MixtureBaseT;
      using typename RpMixtureT::FieldT;

      using MixtureTmplT::polymer;

      // Public member functions

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in) override;

   private:

      /// Use batched FFTs to compute stress? (faster, but 2x memory)
      bool useBatchedFFT_;

      /**
      * Allocate memory for all blocks
      */
      virtual void allocateBlocks() override;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   extern template class MixtureTmpl< Rpg::Polymer<1>, Rpg::Solvent<1> >;
   extern template class MixtureTmpl< Rpg::Polymer<2>, Rpg::Solvent<2> >;
   extern template class MixtureTmpl< Rpg::Polymer<3>, Rpg::Solvent<3> >;
   namespace Rp {
      extern template class Mixture<1, Rpg::Types<1> >;
      extern template class Mixture<2, Rpg::Types<2> >;
      extern template class Mixture<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Mixture<1>;
      extern template class Mixture<2>;
      extern template class Mixture<3>;
   }
}
#endif
