#ifndef RPC_MIXTURE_H
#define RPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/solvers/MixturePrdc.h>     // base class template
#include <rpc/system/Types.h>             // template argument

namespace Pscf {
namespace Rpc {

   // Forward declarations
   template <int D> class Polymer;
   template <int D> class Solvent;

   using namespace Util;
   using namespace Prdc;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from a partial specialization of the template
   * Prdc::MixturePrdc, and has the same public interface as this base
   * class template.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Mixture : public MixturePrdc<D, Polymer<D>, Solvent<D>, Types<D> >
   {

   public:

      /// Direct (parent) base class.
      using MixturePrdcT
         = typename Prdc::MixturePrdc<D, Polymer<D>, Solvent<D>, Types<D> >;

      // Inherited public type name aliases

      using typename MixturePrdcT::MixtureTmplT;
      using typename MixturePrdcT::PolymerT;
      using typename MixturePrdcT::SolventT;
      using typename MixturePrdcT::BlockT;
      using typename MixturePrdcT::PropagatorT;
      using typename MixturePrdcT::FieldT;
      using typename MixturePrdcT::FFTT;
      using typename MixturePrdcT::WaveListT;

      // Inherited public member functions

      using MixturePrdcT::readParameters;
      using MixturePrdcT::associate;
      using MixturePrdcT::allocate;
      using MixturePrdcT::clearUnitCellData;
      using MixturePrdcT::setKuhn;
      using MixturePrdcT::compute;
      using MixturePrdcT::computeStress;
      using MixturePrdcT::hasStress;
      using MixturePrdcT::createBlockCRGrid;

      using MixtureTmplT::polymer;
      using MixtureTmplT::polymerSpecies;
      using MixtureTmplT::solvent;
      using MixtureTmplT::solventSpecies;

      using MixtureBase::nMonomer;
      using MixtureBase::monomer;
      using MixtureBase::nPolymer;
      using MixtureBase::nSolvent;
      using MixtureBase::nBlock;
      using MixtureBase::vMonomer;
      using MixtureBase::isCanonical;

   protected:

      using MixturePrdcT::mesh;
      using MixturePrdcT::ds;

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

   // Explicit instantiation declarations for derived class
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;

} // namespace Rpc
namespace Prdc {

   // Explicit instantiation declarations for base class
   extern template 
   class MixturePrdc<1, Rpc::Polymer<1>, Rpc::Solvent<1>, Rpc::Types<1> >;
   extern template 
   class MixturePrdc<2, Rpc::Polymer<2>, Rpc::Solvent<2>, Rpc::Types<2> >;
   extern template 
   class MixturePrdc<3, Rpc::Polymer<3>, Rpc::Solvent<3>, Rpc::Types<3> >;

} // namespace Prdc
} // namespace Pscf
#endif
