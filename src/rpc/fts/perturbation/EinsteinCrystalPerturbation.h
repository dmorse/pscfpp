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
      using Perturbation<D>::simulator;
      using Perturbation<D>::system;
      using Perturbation<D>::lambda;
      
   protected:
      
      // Inherited protected data members
      using Perturbation<D>::lambda_;
      
   private:
      
      // Spring constant for the einstein crystal. alpha = 1/chi
      double alpha_;
      
      // Reference w field
      DArray< RField<D> > w0_;
      
      // Eigenvector components of the reference w fields
      DArray< RField<D> > wc0_;
      
      // Current Einstein Crystal hamiltonian 
      double ecHamiltonian_;
      
      // Current Block copolymer hamiltonian 
      double bcpHamiltonian_;
      
      // Saved Einstein Crystal hamiltonian  
      double stateEcHamiltonian_;
      
      // Saved Block copolymer hamiltonian 
      double stateBcpHamiltonian_;
      
      // Reference FieldFileName
      std::string referenceFieldFileName_;
      
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
