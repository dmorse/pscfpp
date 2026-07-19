#ifndef RPG_MIXTURE_H
#define RPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/MixtureBase.h>  // base class template
#include <rpg/system/Types.h>        // base class template argument

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class Polymer;
      template <int D, class T> class Solvent;
      template <int D, class T> class Mixture;
   }
}

namespace Pscf {
namespace Rp {

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
   * \see Rp::MixtureBase
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Mixture<D, Rpg::Types<D> > 
     : public Rp::MixtureBase<D, Rpg::Types<D> >
   {

   public:

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

      /// Base class type aliases
      using RpMixtureT = typename Rp::MixtureBase<D, Rpg::Types<D> >;
      using typename RpMixtureT::MixtureTmplT;
      using typename RpMixtureT::CompositionT;
      using typename RpMixtureT::FieldT;

      using MixtureTmplT::polymer;

   private:

      /// Use batched FFTs to compute stress? (faster, but 2x memory)
      bool useBatchedFFT_;

      /**
      * Allocate memory for all blocks
      */
      virtual void allocateBlocks() override;

   };

} // namespace Rp
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class MixtureTmpl< Rp::Polymer<1, Rpg::Types<1> >, Rp::Solvent<1, Rpg::Types<1> > >;
   extern template 
   class MixtureTmpl< Rp::Polymer<2, Rpg::Types<2> >, Rp::Solvent<2, Rpg::Types<2> > >;
   extern template 
   class MixtureTmpl< Rp::Polymer<3, Rpg::Types<3> >, Rp::Solvent<3, Rpg::Types<3> > >;
   namespace Rp {
      extern template class MixtureBase<1, Rpg::Types<1> >;
      extern template class MixtureBase<2, Rpg::Types<2> >;
      extern template class MixtureBase<3, Rpg::Types<3> >;
      extern template class Mixture<1, Rpg::Types<1> >;
      extern template class Mixture<2, Rpg::Types<2> >;
      extern template class Mixture<3, Rpg::Types<3> >;
   }
}
#endif
