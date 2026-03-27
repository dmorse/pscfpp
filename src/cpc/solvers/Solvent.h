#ifndef CPC_SOLVENT_H
#define CPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <prdc/cpu/CField.h>            // member
#include <complex>

// Forward declaration
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Cpc_Solver_Module
   */
   template <int D>
   class Solvent 
    : public SolventSpecies< std::complex<double> >
   {

   public:

      // Typename aliases

      /// Direct base class, a specialization of Pscf::SolventSpecies.
      using Base = SolventSpecies< std::complex<double> >;

      /// Indirect base, a specialization of Pscf::Species.
      using SpeciesT = SolventSpecies< std::complex<double> >;

      // Member functions

      /**
      * Constructor.
      */
      Solvent();

      /**
      * Destructor.
      */
      ~Solvent();

      /**
      * Create an association with the mesh.
      *
      * \param mesh associated Mesh<D> object
      */
      void associate(Mesh<D> const & mesh);

      /**
      * Allocate memory for concentrationf field.
      */
      void allocate();

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Computes monomer concentration field cField, partition function
      * q, and either the solvent volume fraction phi or solvent chemical
      * potential mu, depending on ensemble. The function takes the
      * chemical potential field wField for the relevant monomer type as
      * its only input argument.
      *
      * \param wField  monomer chemical potential field of relevant type.
      */
      void compute(CField<D> const & wField);

      /**
      * Get the monomer concentration field for this solvent.
      */
      CField<D> const & cField() const;

      // Inherited accessor functions
      using SpeciesT::phi;
      using SpeciesT::mu;
      using SpeciesT::q;
      using SpeciesT::ensemble;
      using Base::monomerId;
      using Base::size;

   protected:

      // Inherited protected functions
      using SpeciesT::setQ;

   private:

      /// Complex concentration field for this solvent.
      CField<D> cField_;

      /// Pointer to associated mesh.
      Mesh<D> const *  meshPtr_;

   };

   // Inline member function

   /*
   * Get complex monomer concentration field for this solvent.
   */
   template <int D>
   inline CField<D> const & Solvent<D>::cField() const
   {  return cField_;  }

   // Explicit instantiation declarations
   extern template class Solvent<1>;
   extern template class Solvent<2>;
   extern template class Solvent<3>;

} // namespace Cpc
} // namespace Pscf
#endif
