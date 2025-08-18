#ifndef R1D_MIXTURE_H
#define R1D_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace R1d
{

   // Forward declaration
   class Domain;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A Mixture contains arrays of Polymer and Solvent objects. Each such
   * object can solve the single-molecule statistical mechanics problem 
   * for an ideal gas of the associated species in a set of specified 
   * chemical potential fields, and thereby compute concentrations for
   * all monomer types and single-molecule partition functions for all
   * molecular species. 
   *
   * A Mixture is associated with a Domain object, which models a 1D
   * spatial domain and a spatial discretization mesh. Knowledge of 
   * the domain and discretization is needed to solve the ideal-gas
   * problem.
   *
   * \ref user_param_mixture_page "Parameter File Format"
   * \ingroup R1d_Solver_Module
   */
   class Mixture : public MixtureTmpl<Polymer, Solvent>
   {

   public:

      // Public typedefs

      /**
      * Field type.
      */
      typedef Propagator::FieldT FieldT;

      // Public member functions

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Destructor.
      */
      ~Mixture();

      /**
      * Read all parameters and initialize.
      *
      * This function reads in a complete description of the chemical
      * structure of all species and the composition of the mixture, as
      * well as the target contour length step size ds.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create an association with the domain and allocate memory.
      * 
      * The function parameter domain must have already been initialized, 
      * by reading its parameters from a file, so that the domain 
      * dimensions and number of grid points are known on entry.
      *
      * \param domain associated Domain object (stores address).
      */
      void setDomain(Domain const & domain);

      /**
      * Reset statistical segment length for one monomer type.
      * 
      * This function resets the kuhn or statistical segment length value
      * for a monomer type, and updates the associcated value in every 
      * block of that monomer type.
      *
      * \param monomerId  monomer type id
      * \param kuhn  new value for the statistical segment length
      */
      void setKuhn(int monomerId, double kuhn);

      /**
      * Compute concentrations.
      *
      * This function calls the compute function of every molecular species
      * (polymers and solvents), species, and then adds the resulting block
      * concentration fields for blocks of each type to compute a total 
      * monomer concentration (or volume fraction) for each monomer type.
      * Upon return, values are set for volume fraction and chemical 
      * potential (mu) members of each species, and the concentration field
      * for each Block and Solvent. The total concentration fields for all 
      * monomer types are returned in the cFields output parameter. 
      *
      * The arrays wFields and cFields must each have capacity nMonomer(), 
      * and contain fields that are indexed by monomer type index. 
      *
      * \param wFields  array of chemical potential fields (input)
      * \param cFields  array of monomer concentration fields (output)
      */
      void 
      compute(DArray<FieldT> const & wFields, DArray<FieldT>& cFields);


      // Inherited public member functions with non-dependent names
      using MixtureTmpl< Polymer, Solvent >::nMonomer;
      using MixtureTmpl< Polymer, Solvent >::nPolymer;
      using MixtureTmpl< Polymer, Solvent >::nSolvent;
      using MixtureTmpl< Polymer, Solvent >::nBlock;
      using MixtureTmpl< Polymer, Solvent >::polymer;
      using MixtureTmpl< Polymer, Solvent >::monomer;
      using MixtureTmpl< Polymer, Solvent >::solvent;
      using MixtureTmpl< Polymer, Solvent >::vMonomer;

   protected:

      // Inherited protected member functions with non-dependent names
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Optimal contour length step size.
      double ds_;

      /// Pointer to associated Domain object.
      Domain const * domainPtr_;

      /// Return associated domain by reference.
      Domain const & domain() const;

   };

   // Inline member function

   /*
   * Get Domain by constant reference (private).
   */
   inline Domain const & Mixture::domain() const
   {   
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

} // namespace R1d
} // namespace Pscf
#endif
