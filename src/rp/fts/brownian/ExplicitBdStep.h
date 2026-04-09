#ifndef RP_EXPLICIT_BD_STEP_H
#define RP_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"                    // base class
#include <util/containers/DArray.h>    // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Explicit Euler-Maruyama Brownian dynamics step.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named ExplicitBdStep,
   * that are defined in Rpc and Rpg namespaces and used in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *    - D : dimension
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_ExplicitBdStep_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class ExplicitBdStep : public T::BdStep
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      ExplicitBdStep(typename T::BdSimulator& simulator);

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true iff the compressor converged, false otherwise
      */
      bool step() override;

   protected:

      using BdStepT = typename T::BdStep;
      using BdStepT::system;
      using BdStepT::simulator;
      using BdStepT::vecRandom;

   private:

      // Local copy of w fields
      DArray< typename T::RField > w_;

      // Change in one component of wc
      typename T::RField dwc_;

      // Normal distributed random numbers
      typename T::RField gaussianField_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

}
}
#endif
