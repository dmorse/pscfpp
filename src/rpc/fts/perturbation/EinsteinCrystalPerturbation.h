#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_H

#include "Perturbation.h"      // base class

namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc{
namespace Cpu{
   template <int D> class RField;
}
}
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   template <int D> class Simulator;

   /**
   * Perturbation for Einstein crystal thermodynamic integration method.
   *
   * \see
   * <ul>
   * <li> \ref rpc_EinsteinCrystalPerturbation_page "Manual Page" </li>
   * <li> \ref psfts_perturb_page "Manual Page" </li>
   * </ul>
   *
   * \ingroup Rpc_Fts_Perturbation_Module
   */
   template <int D>
   class EinsteinCrystalPerturbation : public Perturbation<D>
   {

   public:

      /**
      * Constructor.
      */
      EinsteinCrystalPerturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~EinsteinCrystalPerturbation();

      /**
      * Read parameters from archive.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      */
      virtual void setup();

      /**
      * Compute and return the perturbation to the Hamiltonian.
      *
      * \param unperturbedHamiltonian Hamiltonian in absence of perturbation
      */
      virtual double hamiltonian(double unperturbedHamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      */
      virtual void incrementDc(DArray< RField<D> >& dc);
      
      /**
      * Compute and return derivative of free energy.
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
      
      // Inherited functions
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptionalDArray;
      using Perturbation<D>::simulator;
      using Perturbation<D>::system;
      using Perturbation<D>::lambda;
      
   protected:
      
      // Inherited protected data members
      using Perturbation<D>::lambda_;
      
   private:
      
      // Parameters used in Einstein crystal integration
      DArray<double> epsilon_;
      
      // Reference w field
      DArray< RField<D> > w0_;
      
      // Eigenvector components of the reference w fields
      DArray< RField<D> > wc0_;
      
      // Current Einstein crystal Hamiltonian 
      double ecHamiltonian_;
      
      // Current unperturbed Hamiltonian 
      double unperturbedHamiltonian_;
      
      // Saved Einstein crystal Hamiltonian  
      double stateEcHamiltonian_;
      
      // Saved unperturbed Hamiltonian 
      double stateUnperturbedHamiltonian_;

      // Reference FieldFileName
      std::string referenceFieldFileName_;
      
      // Have epsilon values been set ?
      bool hasEpsilon_;
      
      // Compute eigenvector components of the reference field
      void computeWcReference();
      
   };

   #ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
   // Suppress implicit instantiation
   extern template class EinsteinCrystalPerturbation<1>;
   extern template class EinsteinCrystalPerturbation<2>;
   extern template class EinsteinCrystalPerturbation<3>;
   #endif

}
}
#endif
