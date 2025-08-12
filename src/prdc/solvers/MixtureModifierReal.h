#ifndef PRDC_MIXTURE_MODIFIER_REAL_H
#define PRDC_MIXTURE_MODIFIER_REAL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Prdc {

   /**
   * Modifier for parameters of an associated mixture.
   *
   * Template parameter MT is a mixture class type, which is set to
   * a Mixture<D> class in each of the implementations for real periodic
   * fields.
   *
   * A MixtureModifierReal provides an interface for a limited number of
   * non-const operations that modify the parameters of an associated
   * Mixture object. This interface is used to implement sweep and ramp 
   * operations that need to modify these parameters.
   * 
   * Each PSCF implementation for real periodic fields defines a class 
   * template int \<D\> MixtureModifier that is derived from a template 
   * specialization MixtureModifierReal\< Mixture\<D\> \>. The System\<D\> 
   * class provides public access to its Mixture<D> only through a const 
   * reference, but provides a non-const reference to an associated 
   * MixtureModifier. The MixtureModifier template thus defines a set of 
   * operations that modify parameters of the mixture that may be called 
   * via a non-const reference to a System.
   */
   template <class MT>
   class MixtureModifierReal 
   {

   public:

      /// \name Lifetime
      ///@{

      /**
      * Constructor.
      */
      MixtureModifierReal();

      /**
      * Destructor.
      */
      ~MixtureModifierReal();

      /**
      * Create associations with a Mixture.
      *
      * \param mixture  associated Mixture object
      */
      void associate(MT& mixture);

      ///@}
      /// \name Parameter Modification Functions
      ///@{

      /**
      * Set the statistical segment length for one monomer type.
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
      * Set the volume fraction for a polymer species.
      *
      * \param polymerId  polymer species id
      * \param phi  new value for the species volume fraction
      */
      void setPhiPolymer(int polymerId, double phi);

      /**
      * Set the chemical for a polymer species.
      *
      * \param polymerId  polymer species id
      * \param mu  new value for the species chemical potential
      */
      void setMuPolymer(int polymerId, double mu);

      /**
      * Set the length for a block of a block polymer.
      *
      * \param polymerId  polymer species id
      * \param blockId  block identifier
      * \param length  new value for the block length
      */
      void setBlockLength(int polymerId, int blockId, double length);

      /**
      * Set the volume fraction for a solvent species.
      *
      * \param solventId  solvent species id
      * \param phi  new value for the species volume fraction
      */
      void setPhiSolvent(int solventId, double phi);

      /**
      * Set the chemical for a solvent species.
      *
      * \param solventId  solvent species id
      * \param mu  new value for the species chemical potential
      */
      void setMuSolvent(int solventId, double mu);

      /**
      * Set the size (steric volume) for a solvent species.
      *
      * The size is a dimensionless parameter given by the ratio 
      * size = (solvent molecular volume)/vMonomer, where vMonomer
      * is the monomer reference volume.
      *
      * \param solventId  solvent species id
      * \param size  new value for the species size parameter
      */
      void setSolventSize(int solventId, double size);

      /**
      * Set the monomer reference volume for the mixture.
      *
      * \param vMonomer  new value for the monomer reference volume
      */
      void setVMonomer(double vMonomer);

      ///@}

      /**
      * Clear all data that depends on the unit cell parameters.
      *
      * This function marks all private data that depends on the values of
      * the unit cell parameters as invalid, so that it can be recomputed
      * before it is next needed. This function should be called after any 
      * change in unit cell parameters.
      */
      void clearUnitCellData();

   private:

      /**
      * Get the mixture by reference.
      */
      MT& mixture();

      /// Pointer to the asssociated mixture.
      MT* mixturePtr_;

   };

} // namespace Prdc
} // namespace Pscf
#endif
