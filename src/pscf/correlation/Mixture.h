#ifndef PSCF_CORRELATION_MIXTURE_H
#define PSCF_CORRELATION_MIXTURE_H

/*
* PSCF - Mixture Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/correlation/Polymer.h>
#include <util/containers/DArray.h>

// Forward declarations
namespace Pscf {
   template <typename WT> class Composition;
}

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /**
   * Correlations in a homogeneous ideal gas mixture.
   *
   * \ingroup Pscf_Correlation_Module
   */
   template <typename WT>
   class Mixture
   {

   public:

      /**
      * Default constructor.
      */
      Mixture();

      /**
      * Constructor.
      *
      * \param mixture  assocated Composition object
      */
      Mixture(Composition<WT> const & mixture);

      /**
      * Destructor.
      */
      ~Mixture();

      /**
      * Create an association with a Mixture.
      *
      * \param mixture  associated Composition
      */
      void associate(Composition<WT> const& mixture);

      /**
      * Allocate private data structures, set immutable private data.
      *
      * This function may only be called once, after reading the Mixture
      * section of a parameter file. It requires knowledge of immutable 
      * properties of a mixture, such as the number of monomer types, the 
      * number of polymer species, and the monomer type id for each block
      * of each polymer. 
      */
      void allocate();

      /**
      * Set mutable private data.
      *
      * This function calls Pscf::Correlation::Polymer:setup for all
      * polymer species, using the current statistical segment length
      * values. By doing so, it computes properties that depend on the 
      * following mutable properties of a mixture: monomer statistical 
      * segment lengths, species volume fractions, polymer block lengths, 
      * and solvent species sizes. 
      * 
      * This function must be called after allocate and before any 
      * correlation function computations, and may be called repeatedly.
      */
      void setup();

      /**
      * Compute ideal gas correlation functions for a monomer type pair.
      *
      * On return, each element of array "correlations" contains the
      * correlation function for monomers of types ma and mb for a value
      * of squared wavenumber given by the corresponding element of kSq.
      *
      * \param ma  1st monomer type index (in)
      * \param mb  2nd monomer type index (in)
      * \param kSq  array of squared wavenumber values (in)
      * \param correlations  array of correlation function values (out)
      */
      void computeOmega(int ma, int mb,
                        Array<double> const & kSq,
                        Array<double> & correlations) const;

      /**
      * Compute total ideal gas density correlation functions.
      *
      * On return, each element of the array "correlations" contains
      * the total density-density correlation function at a value of
      * squared wavenumber given by the corresponding element of kSq.
      *
      * \param kSq  array of squared wavenumber values (in)
      * \param correlations  array of correlation function values (out)
      */
      void computeOmegaTotal(Array<double> const & kSq,
                             Array<double> & correlations) const;

      /**
      * Return reference to a polymer species descriptor.
      *
      * \param i  index for polymer species
      */
      Correlation::Polymer<WT> const & polymer(int i) const;

      /**
      * Return reference to a descriptor for the associated mixture.
      */
      Composition<WT> const & mixture() const;

      /**
      * Has this Mixture been previously allocated?
      */
      bool isAllocated() const;

   private:

      /// Array of Correlation::Polymer objects
      DArray< Correlation::Polymer<WT> > polymers_;

      /// Pointer to the associated Composition object.
      Composition<WT> const * mixturePtr_;

   };

   // Get a constituent Correlation::Polymer<WT> object by const reference.
   template <typename WT> inline 
   Correlation::Polymer<WT> const & Mixture<WT>::polymer(int i) const
   {
      UTIL_CHECK(polymers_.isAllocated());
      return polymers_[i];
   }

   // Get  a descriptor for the parent mixture by const reference.
   template <typename WT> inline 
   Composition<WT> const & Mixture<WT>::mixture() const
   {  return *mixturePtr_; }

   // Has this object been allocated?
   template <typename WT> inline 
   bool Mixture<WT>::isAllocated() const
   {  return polymers_.isAllocated(); }

   // Explicit instantiation declaration
   extern template class Mixture<double>;

} // namespace Correlation
} // namespace Pscf
#endif
