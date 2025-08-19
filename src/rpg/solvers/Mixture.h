#ifndef RPG_MIXTURE_H
#define RPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/solvers/MixturePrdc.h>     // base class template
#include "Polymer.h"                      // base class parameter
#include "Solvent.h"                      // base class parameter

namespace Pscf {
namespace Rpg {

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture is derived from a partial specialization of the template
   * Prdc::MixturePrdc, and has the same public interface as this base
   * class. See the documentation of this base class template for details.
   *
   * \ref user_param_mixture_page "Manual Page"
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Mixture : public MixturePrdc<D, Polymer<D>, Solvent<D> >
   {

   public:
 
      /// Direct (parent) base class.
      using MixturePrdcT 
         = typename Prdc::MixturePrdc<D, Polymer<D>, Solvent<D> >;

      // Inherited public type name aliases 
      using typename MixturePrdcT::MixtureTmplT;
      using typename MixturePrdcT::PolymerT;
      using typename MixturePrdcT::SolventT;
      using typename MixturePrdcT::BlockT;
      using typename MixturePrdcT::PropagatorT;
      using typename MixturePrdcT::FieldT;
      using typename MixturePrdcT::FFTT;
      using typename MixturePrdcT::WaveListT;

      // Public member functions
      
      /**
      * Constructor.
      */
      Mixture();

      /**
      * Read all parameters and initialize.
      *
      * \param in  input parameter stream 
      */
      void readParameters(std::istream& in) override;

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

      // Private member data
     
      /// Use batched FFTs to compute stress? (faster, but doubles memory)
      bool useBatchedFFT_;

      // Private member functions
      
      /**
      * Set all elements of a field to a common scalar: A[i] = s.
      *
      * \param A  field (LHS)
      * \param s  scalar (RHS)
      */
      virtual void eqS(FieldT& A, double s) const override;

      /**
      * Compound addition assignment for fields : A[i] += B[i].
      *
      * \param A  field (LHS)
      * \param B  field (RHS)
      */
      virtual void addEqV(FieldT& A, FieldT const & B) const override;

      /**
      * Allocate memory for all blocks
      */
      virtual void allocateBlocks() override;  

   };

   // Suppress implicit instantiation
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;

} // namespace Rpg
namespace Prdc {
   // Suppress implicit instantiation of base class 
   extern template class MixturePrdc<1, Rpg::Polymer<1>, Rpg::Solvent<1> >;
   extern template class MixturePrdc<2, Rpg::Polymer<2>, Rpg::Solvent<2> >;
   extern template class MixturePrdc<3, Rpg::Polymer<3>, Rpg::Solvent<3> >;
} // namespace Prdc
} // namespace Pscf
#endif
