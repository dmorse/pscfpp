#ifndef RP_BD_MOVE_H
#define RP_BD_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>  // base class
#include <prdc/field/RField.h>         // member
#include <util/containers/DArray.h>    // member

#include <pscf/backends/TmplDeclare.h>
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Brownian dynamics Monte-Carlo move.
   *
   * A BdMove is a Monte-Carlo move that simply runs a short 
   * Brownian dynamics (BD) simulation, and always accepts the result.
   * This implementation uses the Leimkuhler-Matthews BD step algorithm.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named BdMove, that are
   * defined in Rpc and Rpg namespaces for use in the pscf_rpc pscf_rpg 
   * programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \see \ref rp_BdMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class BdMove : public McMove<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      BdMove(McSimulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~BdMove() = default;


      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Generate a short BD simulation.
      */
      bool move() override;

      /**
      * Do dc derivative components need to be saved before a step?
      *
      * The default implementation returns true.
      *
      * \return true to save, or false otherwise
      */
      bool needsDc() override
      {  return true; }

   protected:

      /**
      * Setup before simulation loop.
      */
      void bdSetup();

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true if compressor converged, false otherwise
      */
      bool bdStep();

      using McMoveT = McMove<D,T>;

      // Inherited protected member functions
      using McMoveT::system;
      using McMoveT::simulator;
      using McMoveT::random;
      using McMoveT::vecRandom;

   private:

      using RFieldT = RField<D,T>;

      // Private data members

      // New field values
      DArray< RFieldT > w_;

      // Random displacements (A)
      DArray< RFieldT > etaA_;

      // Random displacements (B)
      DArray< RFieldT > etaB_;

      // Change in one field component
      RFieldT dwc_;

      // Pointer to new random displacements
      DArray< RFieldT >* etaNewPtr_;

      // Pointer to old random displacements
      DArray< RFieldT >* etaOldPtr_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

      // Number of BD steps per MC move.
      double nStep_;

      // Private member functions

      RFieldT& etaNew(int i)
      {  return (*etaNewPtr_)[i]; }

      RFieldT& etaOld(int i)
      {  return (*etaOldPtr_)[i]; }

      /// Generate new values for etaNew
      void generateEtaNew();

      /// Exchange pointer values for etaNew and etaOld.
      void exchangeOldNew();

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(BdMove)

}
}
#endif
