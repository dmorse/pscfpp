#ifndef RPC_MIXTURE_H
#define RPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/solvers/MixtureReal.h>     // base class template
#include "Polymer.h"                      // base class parameter
#include "Solvent.h"                      // base class parameter

namespace Pscf {
namespace Rpc {

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from a partial specialization of the template
   * Prdc::MixtureReal, and has the same public interface as this base
   * class template.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Mixture : public MixtureReal<D, Polymer<D>, Solvent<D> >
   {

   public:
 
      /// Direct (parent) base class.
      using MixtureRealT 
         = typename Prdc::MixtureReal<D, Polymer<D>, Solvent<D> >;

      // Inherited public type name aliases 
      using typename MixtureRealT::MixtureTmplT;
      using typename MixtureRealT::PolymerT;
      using typename MixtureRealT::SolventT;
      using typename MixtureRealT::BlockT;
      using typename MixtureRealT::PropagatorT;
      using typename MixtureRealT::FieldT;
      using typename MixtureRealT::FFTT;
      using typename MixtureRealT::WaveListT;

      // Inherited public member functions

      using MixtureRealT::readParameters;
      using MixtureRealT::associate;
      using MixtureRealT::allocate;
      using MixtureRealT::clearUnitCellData;
      using MixtureRealT::setKuhn;
      using MixtureRealT::compute;
      using MixtureRealT::computeStress;
      using MixtureRealT::hasStress;
      using MixtureRealT::createBlockCRGrid;

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

      using MixtureRealT::mesh;
      using MixtureRealT::ds;

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

   // Suppress implicit instantiation
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;

} // namespace Rpc
namespace Prdc {
   // Suppress implicit instantiation of base class 
   extern template class MixtureReal<1, Rpc::Polymer<1>, Rpc::Solvent<1> >;
   extern template class MixtureReal<2, Rpc::Polymer<2>, Rpc::Solvent<2> >;
   extern template class MixtureReal<3, Rpc::Polymer<3>, Rpc::Solvent<3> >;
} // namespace Prdc
} // namespace Pscf
#endif
