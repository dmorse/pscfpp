#ifndef PSCF_FH_MIXTURE_H
#define PSCF_FH_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class

#include <pscf/chem/Monomer.h>              // member 
#include <pscf/floryHuggins/FhMolecule.h>   // member 
#include <util/containers/DArray.h>         // member 
#include <util/containers/DMatrix.h>        // member 

// Forward declarations
namespace Pscf {
   class FhInteraction;
   class LuSolver;
   template <typename WT> class Composition;
}

namespace Pscf {

   using namespace Util;

   /**
   * A spatially homogeneous mixture.
   *
   * \ingroup Pscf_FloryHuggins_Module
   */
   class FhMixture : public ParamComposite
   {
   public:

      /**
      * Constructor.
      */
      FhMixture();

      /**
      * Destructor.
      */
      ~FhMixture();

      /// \name Initialization.
      //@{

      /**
      * Read parameters from file and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Initialize to properties of a Composition descriptor.
      *
      * \param mixture  descriptor for a SCFT or FTS mixture
      */
      template <typename WT>
      void initialize(Composition<WT> const& mixture);

      /**
      * Set the number of molecular species and allocate memory.
      *
      * \param nMolecule number of molecular species (polymer and solvent)
      */
      void setNMolecule(int nMolecule);
      
      /**
      * Set the number of monomer types.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      //@}
      /// \name Thermodynamics Computations
      //@{

      /**
      * Set system composition.
      *
      * \param phi array of molecular volume fractions.
      */
      void setComposition(DArray<double> const & phi);

      /**
      * Compute chemical potential from preset composition.
      *
      * Precondition: setComposition must be called prior.
      * Postcondition: Upon return, mu array is set.
      *
      * \param interaction  excess free energy model (input)
      * \param xi  Lagrange multiplier field (input)
      */
      void computeMu(FhInteraction const & interaction, double xi = 0.0);

      /**
      * Compute composition from chemical potentials.
      *
      * \param interaction  excess free energy model (input)
      * \param mu  target molecular chemical potentials (input)
      * \param phi guess of molecular volume fractions (input)
      * \param xi  Lagrange multiplier field (input/output)
      */
      void computePhi(FhInteraction const & interaction, 
                      DArray<double> const & mu, 
                      DArray<double> const & phi, 
                      double&  xi);

      /**
      * Compute Helmholtz free energy and pressure.
      * 
      * Preconditions and postconditions:
      * 
      * \pre  setComposition must be called prior.
      * \pre  computeMu must be called prior.
      * \post fHelmholtz and pressure are set.
      *
      * \param interaction  excess free energy model (input)
      */
      void computeFreeEnergy(FhInteraction const & interaction);

      //@}
      /// \name Accessors 
      //@{
 
      /**
      * Get a molecule object (non-const reference).
      *
      * \param id integer molecule species index (0 <= id < nMolecule)
      */
      FhMolecule& molecule(int id);

      /** 
      * Return chemical potential for one species.
      *
      * \param id integer molecule species index (0 <= id < nMolecule)
      */
      double mu(int id) const;

      /** 
      * Return molecular volume fraction for one species.
      *
      * \param id integer molecule species index (0 <= id < nMolecule)
      */
      double phi(int id) const;

      /** 
      * Return monomer volume fraction for one monomer type.
      *
      * \param id monomer type index (0 <= id < nMonomer)
      */
      double c(int id) const;

      /**
      * Return Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz() const;

      /**
      * Return pressure in units of kT / monomer volume.
      */
      double pressure() const;

      /**
      * Get number of molecule species (polymer + solvent).
      */
      int nMolecule() const;

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      //@}

      /**
      * Validate all data structures.
      *
      * Throw an exception if an error is found.
      */
      void validate() const;

   private:

      /**
      * Array of molecule species solver objects.
      *
      * Array capacity = nMolecule.
      */
      DArray<FhMolecule> molecules_;

      /**
      * Array of molecular chemical potentials. 
      */
      DArray<double> mu_;

      /**
      * Array of molecular volume fractions.
      */
      DArray<double> phi_;

      /**
      * Array of monomer volume fractions.
      */
      DArray<double> c_;

      /**
      * Array of monomer excess chemical potentials.
      */
      DArray<double> w_;

      /**
      * Residual array for used by computePhi function.
      */
      DArray<double> residual_;

      /**
      * Change in input variables (phi, xi)
      */
      DArray<double> dX_;

      /**
      * Derivatives of W with respect to monomer fractions.
      */
      DMatrix<double> dWdC_;

      /**
      * Derivatives of W with respect to molecule fractions.
      */
      DMatrix<double> dWdPhi_;

      /**
      * Jacobian matrix for use by computePhi function.
      */
      DMatrix<double> jacobian_;

      /**
      * Array of old molecular volume fractions, for use in computePhi.
      */
      DArray<double> phiOld_;

      /**
      * Free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Pressure x monomer volume / kT.
      */
      double pressure_;

      /**
      * Pointer to LUSolver.
      */
      LuSolver* solverPtr_;

      /**
      * Number of molecule species (polymers and solvent).
      */
      int nMolecule_;

      /**
      * Number of monomer types (maximum monomer id + 1).
      */
      int nMonomer_;

      /**
      * Initialized false, set true by setComposition().
      */
      bool hasComposition_;

      /**
      * Compute monomer concentrations from phi_.
      */
      void computeC();

      /**
      * Adjust xi to minimize mean-squared residual.
      */
      void adjustXi(DArray<double> const & mu, double& xi); 

      /**
      * Compute residual array and return max error.
      */
      void computeResidual(DArray<double> const & mu,
                           double& error);

   };

   // Inline member functions

   inline FhMolecule& FhMixture::molecule(int id)
   {  
      UTIL_ASSERT(id >= 0);  
      UTIL_ASSERT(id < nMolecule_);  
      return molecules_[id]; 
   }

   inline double FhMixture::mu(int id) const
   {  
      UTIL_ASSERT(id >= 0);  
      UTIL_ASSERT(id < nMolecule_);  
      return mu_[id]; 
   }

   inline double FhMixture::phi(int id) const
   {
      UTIL_ASSERT(id >= 0);  
      UTIL_ASSERT(id < nMolecule_);  
      return phi_[id]; 
   }

   inline double FhMixture::c(int id) const
   {  
      UTIL_ASSERT(id >= 0);  
      UTIL_ASSERT(id < nMonomer_);  
      return c_[id]; 
   }

   inline double FhMixture::fHelmholtz() const
   {  return fHelmholtz_; }

   inline double FhMixture::pressure() const
   {  return pressure_; }

   inline int FhMixture::nMolecule() const
   {  return nMolecule_; }

   inline int FhMixture::nMonomer() const
   {  return nMonomer_; }

   // Explicit instantiation declaration for member function template
   extern template 
   void FhMixture::initialize<double>(Composition<double> const &);

} // namespace Pscf
#endif
