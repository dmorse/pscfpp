#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_H
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_H

#include <rp/fts/perturbation/EinsteinCrystalPerturbation.h> // base class
#include <rpc/system/Types.h>                                // base argument
#include <rpc/fts/perturbation/Perturbation.h>               // indirect base
#include <prdc/field/cpu/RField.h>                           // base member

#if 0
namespace Pscf {
namespace Rpc {

   /**
   * Perturbation for Einstein crystal thermodynamic integration.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the template Rp::EinsteinCrystalPerturbation, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see \ref rp_EinsteinCrystalPerturbation_page "Einstein Crystal"
   * \see \ref psfts_perturb_page "Perturbations"
   * \ingroup Rpc_Fts_Perturbation_Module
   */
   template <int D>
   class EinsteinCrystalPerturbation 
     : public Rp::EinsteinCrystalPerturbation< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      EinsteinCrystalPerturbation(Rp::Simulator<D, Rpc::Types<D> >& simulator);

   };

}
}
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class EinsteinCrystalPerturbation<1, Rpc::Types<1> >;
      extern template 
      class EinsteinCrystalPerturbation<2, Rpc::Types<2> >;
      extern template 
      class EinsteinCrystalPerturbation<3, Rpc::Types<3> >;
   }
}
#endif
