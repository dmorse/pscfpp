#ifndef PSPC_PERTURBATION_TPP
#define PSPC_PERTURBATION_TPP

#include "Perturbation.h"
#include <pspc/simulate/Simulator.h>
#include <prdc/cpu/RField.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   Perturbation<D>::Perturbation(Simulator<D>& simulator)
    : ParamComposite(),
      simulatorPtr_(&simulator)
   {}
   
   /* 
   * Destructor.
   */
   template <int D>
   Perturbation<D>::~Perturbation()
   {}
   
   /*
   * Read parameters from stream, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::readParameters(std::istream& in)
   {}

   /*
   * Setup before simulation, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::setup()
   {}

   /*
   * Compute and return perturbation, default implementation.
   */
   template <int D>
   double Perturbation<D>::hamiltonian()
   { return 0.0; }

   /*
   * Modify functional derivatives, empty default implementation.
   */
   template <int D>
   void Perturbation<D>::incrementDc(DArray< RField<D> > & dc)
   {}

}
}
#endif 
