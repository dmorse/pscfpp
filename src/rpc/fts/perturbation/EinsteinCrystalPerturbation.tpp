#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/System.h>
#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::EinsteinCrystalPerturbation(
                                                   Simulator<D>& simulator)
    : Perturbation<D>(simulator),
      ecHamiltonian_(0.0),
      unperturbedHamiltonian_(0.0),
      stateEcHamiltonian_(0.0),
      stateUnperturbedHamiltonian_(0.0)
   {  setClassName("EinsteinCrystal"); }

   /*
   * Destructor.
   */
   template <int D>
   EinsteinCrystalPerturbation<D>::~EinsteinCrystalPerturbation()
   {}

   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::readParameters(std::istream& in)
   {
      read(in, "lambda", lambda_);
      
      // Optionally read the parameters used in Einstein crystal integration
      const int nMonomer = system().mixture().nMonomer();
      epsilon_.allocate(nMonomer - 1);
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         epsilon_[i] = 0.0;
      }
      readOptionalDArray(in, "epsilon", epsilon_, nMonomer -1);
      
      read(in, "referenceFieldFileName", referenceFieldFileName_);
   }

   /*
   * Setup before simulation.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::setup()
   {
      const int nMonomer = system().mixture().nMonomer();
      const IntVec<D>
      meshDimensions = system().domain().mesh().dimensions();

      UTIL_CHECK(nMonomer == 2);

      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(meshDimensions);
         wc0_[i].allocate(meshDimensions);
      }
      
      /* 
      * If the user does not input values for epsilon, set 
      * the values of epsilon_ to the nontrivial -1.0 * eigenvalues 
      * of the projected chi matrix by default.
      */ 
      if (epsilon_[0] == 0){
         for (int i = 0; i < nMonomer - 1 ; ++i) {
            epsilon_[i] = -1.0 * simulator().chiEval(i);
         }
      }

      // Read in reference field from a file
      UnitCell<D> tempUnitCell;
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGrid(referenceFieldFileName_,
                              w0_, tempUnitCell);

      // Compute eigenvector components of the reference field
      computeWcReference();
   }

   /*
   * Compute and return perturbation to Hamiltonian.
   */
   template <int D>
   double
   EinsteinCrystalPerturbation<D>::hamiltonian(double unperturbedHamiltonian)
   {
      // Compute Einstein crystal Hamiltonian
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      double prefactor, w;
      ecHamiltonian_ = 0.0;

      for (int j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & Wc = simulator().wc(j);
         prefactor = double(nMonomer)/(2 * epsilon_[j] );
         for (int i = 0; i < meshSize; ++i) {
            w = Wc[i] - wc0_[j][i];
            ecHamiltonian_ += prefactor*w*w;
         }
      }

      // Normalize EC hamiltonian to equal a value per monomer
      ecHamiltonian_ /= double(meshSize);

      // Compute EC hamiltonian of system
      ecHamiltonian_ *= nMonomerSystem;

      // Obtain block copolymer hamiltonian
      unperturbedHamiltonian_ = unperturbedHamiltonian;

      return (1.0 - lambda_)*(ecHamiltonian_ - unperturbedHamiltonian_);
   }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void
   EinsteinCrystalPerturbation<D>::incrementDc(DArray< RField<D> > & dc)
   {
      const int meshSize = system().domain().mesh().size();
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      double DcEC, prefactor;

      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         RField<D>& Dc = dc[i];
         RField<D> const & Wc = simulator().wc(i);

         // Loop over grid points
         for (int k = 0; k < meshSize; ++k) {
            prefactor = double(nMonomer) / (epsilon_[i] * vMonomer);
            DcEC = prefactor * (Wc[k] - wc0_[i][k]);

            // Compute composite derivative
            Dc[k] += (1.0 - lambda_) * (DcEC - Dc[k]);
         }
      }
   }

   /*
   * Compute and return derivative of free energy with respect to lambda.
   */
   template <int D>
   double EinsteinCrystalPerturbation<D>::df(){
      return unperturbedHamiltonian_ - ecHamiltonian_;
   }

   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::saveState(){
      stateEcHamiltonian_ = ecHamiltonian_;
      stateUnperturbedHamiltonian_ = unperturbedHamiltonian_;
   }

   /*
   * Save any required internal state variables.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::restoreState(){
      ecHamiltonian_ = stateEcHamiltonian_;
      unperturbedHamiltonian_ = stateUnperturbedHamiltonian_;
   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs of the projected chi matrix as a basis.
   */
   template <int D>
   void EinsteinCrystalPerturbation<D>::computeWcReference()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Loop over eigenvectors (j is an eigenvector index)
      for (j = 0; j < nMonomer; ++j) {

         // Loop over grid points to zero out field wc_[j]
         RField<D>& Wc = wc0_[j];
         for (i = 0; i < meshSize; ++i) {
            Wc[i] = 0.0;
         }

         // Loop over monomer types (k is a monomer index)
         for (k = 0; k < nMonomer; ++k) {
            double vec = simulator().chiEvecs(j, k)/double(nMonomer);
            
            // Loop over grid points
            RField<D> const & Wr = w0_[k];
            for (i = 0; i < meshSize; ++i) {
               Wc[i] += vec*Wr[i];
            }
         }
      }

      #if 0
      Log::file() << "wc " << wc0_.capacity() << "\n";
      for (i = 0; i < 10; ++i) {
         Log::file() << "wc_1 " << wc0_[0][i] << "\n";
         Log::file() << "wc_2 " << wc0_[1][i] << "\n";
      }
      #endif
   }

}
}
#endif
