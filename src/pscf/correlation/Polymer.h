#ifndef PSCF_CORRELATION_POLYMER_H
#define PSCF_CORRELATION_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header includes
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/GArray.h>

// Forward declarations
namespace Pscf {
   template <typename WT> class PolymerSpecies;
}

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /**
   * Intramolecular correlation analysis for one polymer Species.
   *
   * \ingroup Pscf_Correlation_Module
   */
   template <typename WT>
   class Polymer
   {

   public:

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Constructor (creates association with a PolymerSpecies).
      *
      * \param polymer  associated PolymerSpecies object
      */
      Polymer(PolymerSpecies<WT> const & polymer);

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Create association with a PolymerSpecies.
      *
      * \param polymer  associated PolymerSpecies object
      */
      void associate(PolymerSpecies<WT> const & polymer);

      /**
      * Allocate memory and initialize immutable data.
      *
      * This function must be called exactly once, at any time after the
      * Mixture block of the parameter file has been processes. It performs
      * memory allocation and initialization operations that require
      * knowlege of the number of monomer types in the mixture, as well
      * as the number of blocks and the monomer type index of each block
      * in the associated polymer species.
      *
      * \pre An association with a PolymerSpecies must exist.
      *
      * \param nMonomer  number of monomer types in the mixture
      */
      void allocate(int nMonomer);

      /**
      * Set private mutable data.
      *
      * This function must be called after allocate and before performing
      * any computations of correlation functions. It sets internal
      * variables that depend on the polymer volume fraction phi, block
      * lengths and statistical segment lengths. These variables are all
      * mutable insofar as they may all change during a SCFT sweep or FTS
      * ramp, and phi may change in an open statistical ensemble. This
      * function may be called multiple times.
      *
      * \pre The allocate function must have been called previously.
      *
      * \param kuhn  array of kuhn lengths, indexed by monomer type id
      */
      void setup(Array<double> const & kuhn);

      /**
      * Compute intramolecular correlation function for a pair of blocks.
      *
      * This function evaluates and returns the intramolecular correlation
      * function Omega(k) for a pair of blocks with block indices ia and
      * ib at a value of the squared wavenumber given by parameter "kSq".
      * Monomer type indexes ia and ib may be either equal or unequal.
      * The return value is given by the product of a dimensionless
      * single-chain correlation function omega(k) and the input parameter
      * "prefactor". To obtain the contribution of this function to a
      * correlation function for a homogeneous ideal gas mixture, set
      * the prefactor parameter equal to the species chain concentration,
      * prefactor = phi/(v*totalLength), where v is the monomer reference
      * volume.
      *
      * \pre The setup member function must have been called previously.
      *
      * \pre Array parameters kSq and correlations must be allocated with
      * equal nonzero capacities.
      *
      * \param ia  block index of first block
      * \param ib  block index of second block
      * \param prefactor  prefactor multiplying omega(k)
      * \param kSq  squared wavenumber value
      * \return calculated value of Omega_{ia,ib}(k)
      */
      double computeOmega(int ia, int ib, double prefactor, double kSq)
      const;

      /**
      * Compute intramolecular correlation function for a pair of blocks.
      *
      * This function evaluates the intramolecular correlation function
      * Omega for a pair of blocks with block indices ia and ib at a list
      * of values for the squared wavenumber given by values of elements
      * of input array parameter "kSq". Results are added to corresponding
      * elements of array "correlations". The value for each wavenumber is
      * given by the product of a dimensionless single-chain correlation
      * function omega(k) and the input parameter "prefactor". To obtain
      * the contribution of this species to a correlation function for an
      * ideal gas mixture, set prefactor = phi/(v*totalLength), where v
      * is the monomer reference volume.
      *
      * Resulting values of Omega are added to prior values of elements of
      * the array correlation. This is done to simplify calculations in
      * which results of multiple polymer species are added to obtain a
      * correlation function for a mixture. All elements of the array
      * correlations should thus be set to zero at the beginning any such
      * calculation.
      *
      * \pre The setup member function must have been called previously.
      *
      * \pre Array parameters kSq and correlations must be allocated with
      * equal nonzero capacities.
      *
      * \param ia  block index of first block
      * \param ib  block index of second block
      * \param prefactor  prefactor multiplying omega(q)
      * \param kSq  array of squared wavenumbers (in)
      * \param correlation  array of correlation function values (out)
      */
      void computeOmega(int ia, int ib, double prefactor,
                        Array<double> const & kSq,
                        Array<double> & correlation) const;

      /**
      * Compute total intramolecular correlation function.
      *
      * This function computes the a total density-density correlation
      * function for a list of values for the squared wavenumber given
      * as elements of input array parameter "kSq". Results are added
      * to corresponding elements of array parameter "correlations".
      * On entry, the arrays kSq and correlations must be allocated and
      * equal capacities.
      *
      * Each computed value is given by the product of a dimensionless
      * single-chain correlation function omega(k) and the "prefactor"
      * input parameter. To obtain the contribution of this species to
      * the total density-density correlation function for an ideal gas
      * mixture, set prefactor = phi/(v*totalLength), where v is the
      * monomer reference volume.  The zero wavenumber limiting value
      * is given by prefactor*totalLength*totalLength.
      *
      * Resulting values of Omega are added to elements of output array
      * correlation, as for the function computeOmega.
      *
      * \pre The setup member function must have been called previously.
      *
      * \pre Array parameters kSq and correlations must be allocated with
      * equal nonzero capacities.
      *
      * \param prefactor  prefactor of single-chain correlation functions
      * \param kSq  array of squared wavenumbers (in)
      * \param correlations  array of correlation function values (out)
      */
      void computeOmegaTotal(double prefactor,
                             Array<double> const & kSq,
                             Array<double> & correlations) const;

      /**
      * Return the volume fraction of this polymer species.
      */
      double phi() const;

      /**
      * Return the sum of lengths of all blocks in this polymer species.
      */
      double totalLength() const;

      /**
      * Get the number of blocks in the associated polymer species.
      */
      int nBlock() const;

      /**
      * Get a GArray of block ids for all blocks of one monomer type.
      *
      * This function returns a GArray<int> container that contains a
      * list of block indices for all blocks of the specified monomer
      * type i.  The number of blocks of a monomer type i is given by
      * the return value of the "size()" member function of this
      * container.  If this polymer species does not contain any blocks
      * of monomer type i, this size will be zero.
      *
      * \param i  monomer type index
      */
      GArray<int> const & blockIds(int i) const;

      /**
      * Get the length of each block.
      *
      * This function returns the length of block i for the thread model,
      * or the floating point representation of nBead for the the bead
      * model.
      *
      * \param i  block index
      */
      double length(int i) const;

      /**
      * Get the mean-squared length of path between blocks i and j.
      *
      * This function returns the sum of the mean-squared end-to-end
      * lengths of all of the blocks (for the thread model) or bonds
      * (for the bead model) along the path connecting the two blocks
      * with block indices i and j. For the bead model, this includes
      * the two "half bonds" that connect adjacent bonds along this
      * path to each shared vertex.
      *
      * \param i  block index of 1st block.
      * \param j  block index of 2nd block.
      */
      double rSq(int i, int j) const;

   private:

      /*
      * Block lengths, indexed by block index.
      *
      * For the thread model, each element contains the length parameter
      * for associated block. For the bead model, each element contains
      * a double precision representation of nBead, the number of beads.
      */
      DArray<double> length_;

      /*
      * Monomer statistical segment lengths, indexed by block index.
      */
      DArray<double> kuhn_;

      /*
      * Symmetric matrix of values of rSq for paths between blocks.
      *
      * Diagonal elements are zero.
      */
      DMatrix<double> rSq_;

      /*
      * Identifiers for blocks of specified monomer type.
      *
      * Elements of the outer DArray are indexed by monomer type id.
      * Each element is a GArray<int> container containing a list of
      * blocks ids for all blocks of the relevant monomer type, if any.
      */
      DArray< GArray<int> > blockIds_;

      /*
      * Pointer to the associated PolymerSpecies object.
      */
      PolymerSpecies<WT> const * speciesPtr_;

      /*
      * Volume fraction of this polymer.
      */
      double phi_;

      /*
      * Sum of lengths of all blocks in the polymer.
      */
      double totalLength_;

      /*
      * Number of blocks in this polymer.
      */
      int nBlock_;

      /*
      * Number of monomer types in the mixture.
      */
      int nMonomer_;

      /**
      * Return reference to associated polymer species.
      */
      PolymerSpecies<WT> const & species() const;

   };

   // Public inline member functions

   // Get the volume fraction of this polymer species.
   template <typename WT> inline
   double Polymer<WT>::phi() const
   {  return phi_; }

   // Get the sum of the lengths of all blocks in this polymer species.
   template <typename WT> inline
   double Polymer<WT>::totalLength() const
   {  return totalLength_; }

   // Get the number of blocks in this polymer species.
   template <typename WT> inline
   int Polymer<WT>::nBlock() const
   {  return nBlock_; }

   // Get the list of blocks ids for a specified monomer type.
   template <typename WT> inline
   GArray<int> const & Polymer<WT>::blockIds(int i) const
   {  return blockIds_[i]; }

   // Get the length of a specified block.
   template <typename WT> inline
   double Polymer<WT>::length(int i) const
   {  return length_[i]; }

   // Get the mean-squared length of the path connecting two blocks.
   template <typename WT> inline
   double Polymer<WT>::rSq(int i, int j) const
   {  return rSq_(i, j); }

   // Private inline member function

   // Get the associated PolymerSpecies molecule descriptor object.
   template <typename WT> inline
   PolymerSpecies<WT> const & Polymer<WT>::species() const
   {  return *speciesPtr_; }

   // Explicit specialization declaration
   extern template class Polymer<double>;

} // namespace Correlation
} // namespace Pscf
#endif
