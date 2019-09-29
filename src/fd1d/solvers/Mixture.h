#ifndef FD1D_MIXTURE_H
#define FD1D_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Fd1d
{

   class Domain;

   /**
   * Mixture of polymers and solvents.
   *
   * A Mixture contains a list of Polymer and Solvent objects. Each
   * such object can solve the single-molecule statistical mechanics 
   * problem for an ideal gas of the associated species in a set of
   * specified chemical potential fields, and thereby compute 
   * concentrations and single-molecule partition functions. A
   * Mixture is thus both a chemistry descriptor and an ideal-gas 
   * solver.
   *
   * A Mixture is associated with a Domain object, which models a
   * spatial domain and a spatial discretization. Knowledge of the
   * domain and discretization is needed to solve the ideal-gas
   * problem.
   *
   * \ingroup Fd1d_Solver_Module
   */
   class Mixture : public MixtureTmpl<Polymer, Solvent>
   {

   public:

      // Public typedefs

      /**
      * Monomer chemical potential field type.
      */
      typedef Propagator::WField WField;

      /**
      * Monomer concentration or volume fraction field type.
      */
      typedef Propagator::CField CField;

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
      * This function reads in a complete description of
      * the chemical composition and structure of all species,
      * as well as the target contour length step size ds.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create an association with the domain and allocate memory.
      * 
      * The domain parameter must have already been initialized, 
      * e.g., by reading its parameters from a file, so that the
      * domain dimensions are known on entry.
      *
      * \param domain associated Domain object (stores address).
      */
      void setDomain(Domain const & domain);

      /**
      * Compute concentrations.
      *
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration
      * fields for blocks of each type to compute a total monomer
      * concentration (or volume fraction) for each monomer type.
      * Upon return, values are set for volume fraction and chemical 
      * potential (mu) members of each species, and for the 
      * concentration fields for each Block and Solvent. The total
      * concentration for each monomer type is returned in the
      * cFields output parameter. 
      *
      * The arrays wFields and cFields must each have size nMonomer(),
      * and contain fields that are indexed by monomer type index. 
      *
      * \param wFields array of chemical potential fields (input)
      * \param cFields array of monomer concentration fields (output)
      */
      void 
      compute(DArray<WField> const & wFields, DArray<CField>& cFields);

      /**
      * Get monomer reference volume.
      */
      double vMonomer() const;

   private:

      /// Monomer reference volume (set to 1.0 by default).
      double vMonomer_;

      /// Optimal contour length step size.
      double ds_;

      /// Pointer to associated Domain object.
      Domain const * domainPtr_;

      /// Return associated domain by reference.
      Domain const & domain() const;

   };

   // Inline member function

   /*
   * Get monomer reference volume (public).
   */
   inline double Mixture::vMonomer() const
   {  return vMonomer_; }

   /*
   * Get Domain by constant reference (private).
   */
   inline Domain const & Mixture::domain() const
   {   
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
