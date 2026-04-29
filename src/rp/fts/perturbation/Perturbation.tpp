#ifndef RP_PERTURBATION_TPP
#define RP_PERTURBATION_TPP

//Headers to be included in derived class implementation
// T::Simulator
// T::RField

#include "Perturbation.h"
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   Perturbation<D,T>::Perturbation(typename T::Simulator& simulator)
    : ParamComposite(),
      lambda_(1.0),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {}

   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D, class T>
   void Perturbation<D,T>::readParameters(std::istream& in)
   {}

   /*
   * Setup before simulation, empty default implementation.
   */
   template <int D, class T>
   void Perturbation<D,T>::setup()
   {}

   /*
   * Compute and return Hamiltonian perturbation, default implementation.
   */
   template <int D, class T>
   double Perturbation<D,T>::hamiltonian(double unperturbedHamiltonian)
   {  return 0.0; }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D, class T>
   void Perturbation<D,T>::incrementDc(DArray< typename T::RField > & dc)
   {}

   /*
   * Save any internal state variables.
   */
   template <int D, class T>
   void Perturbation<D,T>::saveState()
   {}

   /*
   * Restore any saved internal variables.
   */
   template <int D, class T>
   void Perturbation<D,T>::restoreState()
   {}

   /*
   * Compute and return derivative of free energy.
   */
   template <int D, class T>
   double Perturbation<D,T>::df()
   {  return 0.0; }

   /*
   * Set a new value for the lambda_ parameter.
   */
   template <int D, class T>
   void Perturbation<D,T>::setLambda(double lambda)
   {  lambda_ = lambda; }

}
}
#endif
