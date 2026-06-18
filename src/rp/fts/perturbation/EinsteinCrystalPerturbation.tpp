#ifndef RP_EINSTEIN_CRYSTAL_PERTURBATION_TPP
#define RP_EINSTEIN_CRYSTAL_PERTURBATION_TPP

// Derived classes must include headers for:
//  T::Simulator
//  T::System
//  T::Mixture
//  T::Domain
//  T::VecOp
//  T::Reduce

#include "EinsteinCrystalPerturbation.h"
#include <prdc/crystal/UnitCell.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   EinsteinCrystalPerturbation<D,T>::EinsteinCrystalPerturbation(
                                      typename T::Simulator& simulator)
    : PerturbationT(simulator),
      ecHamiltonian_(0.0),
      unperturbedHamiltonian_(0.0),
      stateEcHamiltonian_(0.0),
      stateUnperturbedHamiltonian_(0.0),
      hasEpsilon_(false)
   {  ParamComposite::setClassName("EinsteinCrystalPerturbation"); }

   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D, class T>
   void EinsteinCrystalPerturbation<D,T>::readParameters(std::istream& in)
   {
      ParamComposite::read(in, "lambda", lambda_);

      // Allocate and initialize epsilon_ array
      const int nMonomer = system().mixture().nMonomer();
      epsilon_.allocate(nMonomer - 1);
      for (int i = 0; i < nMonomer - 1 ; ++i) {
         epsilon_[i] = 0.0;
      }

      // Optionally read array of epsilon parameters
      int nc = nMonomer - 1;
      hasEpsilon_ = ParamComposite::readOptionalDArray(in, "epsilon", 
                                                  epsilon_, nc).isActive();

      ParamComposite::read(in, "fieldFileName", fieldFileName_);
   }

   /*
   * Setup before simulation.
   */
   template <int D, class T>
   void EinsteinCrystalPerturbation<D,T>::setup()
   {
      const int nMonomer = system().mixture().nMonomer();
      const IntVec<D> meshDimensions
                      = system().domain().mesh().dimensions();

      // Allocate memory for reference field
      w0_.allocate(nMonomer);
      wc0_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w0_[i].allocate(meshDimensions);
         wc0_[i].allocate(meshDimensions);
      }
      dw_.allocate(meshDimensions);

      /*
      * If the user did not input epsilon_ values, set values to -1.0
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
      UnitCell<D> tempUnitCell;
      FieldIo<D,T> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGrid(fieldFileName_, w0_, tempUnitCell);

      // Compute eigenvector components of the reference field
      computeWcReference();
   }

   /*
   * Compute and return perturbation to Hamiltonian.
   */
   template <int D, class T>
   double
   EinsteinCrystalPerturbation<D,T>::hamiltonian(
                                          double unperturbedHamiltonian)
   {
      // Set unperturbedHamiltonian_ member variable
      unperturbedHamiltonian_ = unperturbedHamiltonian;

      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      Domain<D,T> const & domain = system().domain();
      const int meshSize = domain.mesh().size();
      const double vSystem  = domain.unitCell().volume();
      const double nMonomerSystem = vSystem / vMonomer;

      // Compute Einstein crystal Hamiltonian
      double prefactor;
      ecHamiltonian_ = 0.0;
      for (int j = 0; j < nMonomer - 1; ++j) {
         VecOp::subVV(dw_, simulator().wc(j), wc0_[j]);
         prefactor = double(nMonomer)/(2.0 * epsilon_[j]);
         ecHamiltonian_ += prefactor * Reduce::sumSq(dw_);
      }
      ecHamiltonian_ /= double(meshSize);
      ecHamiltonian_ *= nMonomerSystem;

      return (1.0 - lambda_)*(ecHamiltonian_ - unperturbedHamiltonian_);
   }

   /*
   * Modify functional derivatives.
   */
   template <int D, class T>
   void
   EinsteinCrystalPerturbation<D,T>::incrementDc(DArray<RFieldT> & dc)
   {
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      double prefactor;

      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         prefactor = double(nMonomer) / (epsilon_[i] * vMonomer);

         // dw_ = prefactor * (wc - wc0_)
         VecOp::subVV(dw_, simulator().wc(i), wc0_[i]);
         VecOp::mulEqS(dw_, prefactor);

         // dc = dc * lambda_ + dw_ * (1 - lambda_)
         VecOp::mulEqS(dc[i], lambda_);
         VecOp::addEqVc(dc[i], dw_, 1.0 - lambda_);
      }
   }

   /*
   * Compute and return derivative of free energy with respect to lambda.
   */
   template <int D, class T>
   double EinsteinCrystalPerturbation<D,T>::df()
   {  return unperturbedHamiltonian_ - ecHamiltonian_; }

   /*
   * Save any required internal state variables.
   */
   template <int D, class T>
   void EinsteinCrystalPerturbation<D,T>::saveState()
   {
      stateEcHamiltonian_ = ecHamiltonian_;
      stateUnperturbedHamiltonian_ = unperturbedHamiltonian_;
   }

   /*
   * Save any required internal state variables.
   */
   template <int D, class T>
   void EinsteinCrystalPerturbation<D,T>::restoreState()
   {
      ecHamiltonian_ = stateEcHamiltonian_;
      unperturbedHamiltonian_ = stateUnperturbedHamiltonian_;
   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs of the projected chi matrix as a basis.
   */
   template <int D, class T>
   void EinsteinCrystalPerturbation<D,T>::computeWcReference()
   {
      const int nMonomer = system().mixture().nMonomer();
      int i, j;

      // Loop over eigenvectors
      for (i = 0; i < nMonomer; ++i) {

         // Initialize to zero
         VecOp::eqS(wc0_[i], 0.0);

         // Loop over monomer types
         for (j = 0; j < nMonomer; ++j) {
            double vec = simulator().chiEvecs(i, j)/double(nMonomer);
            VecOp::addEqVc(wc0_[i], w0_[j], vec);
         }
      }
   }

}
}
#endif
