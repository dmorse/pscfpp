#ifndef RP_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RP_EINSTEIN_CRYSTAL_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <prdc/field/RField.h>                 // member
#include <util/containers/DArray.h>            // member

#include <pscf/backends/TmplDeclare.h>
#include <iostream>

// Forward declaration
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Perturbation for Einstein crystal thermodynamic integration.
   *
   * Specializations of this class template are used as base classes for two
   * analogous class templates, also named EinsteinCrystalPerturbation,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \see \ref rp_EinsteinCrystalPerturbation_page "Einstein Crystal"
   * \see \ref psfts_perturb_page "Perturbations"
   * \ingroup Rp_Fts_Perturbation_Module
   */
   template <int D, class T>
   class EinsteinCrystalPerturbation : public Perturbation<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      EinsteinCrystalPerturbation(Simulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~EinsteinCrystalPerturbation() = default;

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in  input parameter file stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      */
      virtual void setup();

      /**
      * Compute and return the perturbation to the Hamiltonian.
      *
      * \param unperturbedHamiltonian  Hamiltonian without perturbation
      */
      virtual double hamiltonian(double unperturbedHamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      *
      * \param dc  functional derivatives of Hamiltonian (in/out)
      */
      virtual void incrementDc(DArray< RField<D,T> >& dc);

      /**
      * Compute and return derivative of free energy w/ respect to lambda.
      */
      virtual double df();

      /**
      * Save any required internal state variables.
      */
      virtual void saveState();

      /**
      * Restore any required internal state variables.
      */
      virtual void restoreState();

   protected:

      // Alias for base class.
      using PerturbationT = Perturbation<D,T>;

      // Inherited protected members (selected).
      using Perturbation<D,T>::lambda_;
      using Perturbation<D,T>::simulator;
      using Perturbation<D,T>::system;

   private:

      // Parameters used in Einstein crystal integration
      DArray<double> epsilon_;

      // Reference w field
      DArray< RField<D,T> > w0_;

      // Eigenvector components of the reference w fields
      DArray< RField<D,T> > wc0_;

      // Work space
      RField<D,T> dw_;

      // Current Einstein crystal Hamiltonian
      double ecHamiltonian_;

      // Current unperturbed Hamiltonian
      double unperturbedHamiltonian_;

      // Saved Einstein crystal Hamiltonian
      double stateEcHamiltonian_;

      // Saved unperturbed Hamiltonian
      double stateUnperturbedHamiltonian_;

      // Reference field file name
      std::string fieldFileName_;

      // Have epsilon values been set?
      bool hasEpsilon_;

      // Compute eigenvector components of the reference field.
      void computeWcReference();

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(EinsteinCrystalPerturbation)

}
}
#endif
