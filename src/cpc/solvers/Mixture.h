#ifndef CPC_MIXTURE_H
#define CPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <cp/solvers/Mixture.h>        // base class template
#include <cpc/system/Types.h>          // base class template argument

namespace Pscf {
namespace Cpc {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from a partial specialization of the template
   * Cp::Mixture, and has the same public interface as this base class
   * template.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Cpc_Solver_Module
   */
   template <int D>
   class Mixture 
    : public Cp::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >
   {

   public:

      /// Base classes
      using CpMixtureT
         = typename Cp::Mixture<D, Polymer<D>, Solvent<D>, Types<D> >;
      using typename CpMixtureT::MixtureTmplT;
      using typename CpMixtureT::MixtureBaseT;

      using typename CpMixtureT::FieldT;

      using MixtureTmplT::polymer;
      using MixtureBaseT::nPolymer;

   protected:

      using CpMixtureT::mesh;
      using CpMixtureT::ds;

   private:

      /**
      * Set all elements of a field to a common scalar: A[i] = s.
      *
      * \param A  field (LHS)
      * \param s  scalar (RHS)
      */
      void eqS(FieldT& A, double s) const override;

      /**
      * Compound addition assignment for fields : A[i] += B[i].
      *
      * \param A  field (LHS)
      * \param B  field (RHS)
      */
      void addEqV(FieldT& A, FieldT const & B) const override;

      /**
      * Allocate memory for all blocks
      */
      void allocateBlocks() override;

   };

} // namespace Cpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Cp {
      extern template 
      class Mixture<1, Cpc::Polymer<1>, Cpc::Solvent<1>, Cpc::Types<1> >;
      extern template 
      class Mixture<2, Cpc::Polymer<2>, Cpc::Solvent<2>, Cpc::Types<2> >;
      extern template 
      class Mixture<3, Cpc::Polymer<3>, Cpc::Solvent<3>, Cpc::Types<3> >;
   }
   namespace Cpc {
      extern template class Mixture<1>;
      extern template class Mixture<2>;
      extern template class Mixture<3>;
   }
}
#endif
