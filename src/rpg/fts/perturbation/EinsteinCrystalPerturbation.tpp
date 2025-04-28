#ifndef RPG_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RPG_EINSTEIN_CRYSTAL_PERTURBATION_TPP

#include "EinsteinCrystalPerturbation.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/VecOpFts.h>
#include <rpg/System.h>
#include <prdc/cuda/resources.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc::Cuda;

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
      stateUnperturbedHamiltonian_(0.0),
      hasEpsilon_(false)
   {  setClassName("EinsteinCrystalPerturbation"); }

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

      // Allocate and initialize epsilon_ array
      const int nMonomer = system().mixture().nMonomer();
      epsilon_.allocate(nMonomer - 1);
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         epsilon_[i] = 0.0;
      }

      // Optionally read the parameters used in Einstein crystal integration
      hasEpsilon_ =
        readOptionalDArray(in, "epsilon", epsilon_, nMonomer-1).isActive();

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

      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(meshDimensions);
         wc0_[i].allocate(meshDimensions);
      }

      /*
      * If the user did not input values for epsilon, set values to -1.0
      * times the nontrivial eigenvalues of the projected chi matrix.
      */
      if (!hasEpsilon_){
         for (int i = 0; i < nMonomer - 1 ; ++i) {
            epsilon_[i] = -1.0 * simulator().chiEval(i);
         }
      }

      // Check that all epsilon values are positive
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         UTIL_CHECK(epsilon_[i] > 0.0);
      }

      // Read in reference field from a file
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGrid(referenceFieldFileName_,
                              w0_, system().domain().unitCell());

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
      const IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      double prefactor;
      RField<D> wcs;
      wcs.allocate(meshDimensions);
      ecHamiltonian_ = 0.0;

      for (int j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & Wc = simulator().wc(j);
         prefactor = double(nMonomer)/(2.0 * epsilon_[j]);
         VecOp::subVV(wcs, Wc, wc0_[j]); // wcs = Wc - wc0_[j]
         double wSqure = 0;
         wSqure = Reduce::innerProduct(wcs, wcs);
         ecHamiltonian_ += prefactor * wSqure;
      }

      // Normalize EC hamiltonian to equal a value per monomer
      ecHamiltonian_ /= double(meshSize);

      // Compute EC hamiltonian of system
      ecHamiltonian_ *= nMonomerSystem;

      // Set unperturbedHamiltonian_ member variable
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
      const IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      double prefactor;
      RField<D> DcEC;
      DcEC.allocate(meshDimensions);

      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         RField<D>& Dc = dc[i];
         RField<D> const & Wc = simulator().wc(i);
         prefactor = double(nMonomer) / (epsilon_[i] * vMonomer);

         // Compute EC derivative
         // DcEC = prefactor * (Wc - wc0_);
         VecOp::subVV(DcEC, Wc, wc0_[i]);
         VecOp::mulEqS(DcEC, prefactor);

         // Compute composite derivative
         // Dc = (Dc * lambda_) + (DcEC * (1-lambda_))
         VecOp::addVcVc(Dc, Dc, lambda_, DcEC, 1 - lambda_);
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
      int i, j;

      // Loop over eigenvectors (i is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {

         // Loop over grid points to zero out field wc0_[i]
         RField<D>& Wc = wc0_[i];

         VecOp::eqS(Wc, 0.0); // initialize to 0

         // Loop over monomer types (j is a monomer index)
         for (j = 0; j < nMonomer; ++j) {
            cudaReal vec;
            vec = (cudaReal) simulator().chiEvecs(i, j)/nMonomer;

            // Loop over grid points
            VecOp::addEqVc(Wc, w0_[j], vec);
         }

      }
   }

}
}
#endif
