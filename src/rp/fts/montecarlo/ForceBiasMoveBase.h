#ifndef RP_FORCE_BIAS_MOVE_BASE_H
#define RP_FORCE_BIAS_MOVE_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>    // base class
#include <util/containers/DArray.h>      // member
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * ForceBiasMoveBase attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler-Maruyama Brownian
   * dynamics step.
   *
   * Because the probability of attempting a move is not equal to that
   * of generating the reverse move, the acceptance criterion used in
   * the move() function must take into account the ratio of generation
   * probabilities.
   *
   * Template parameters:
   *
   *    - D : dimension of space (D=1, 2, or 3)
   *    - T : Types class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class ForceBiasMoveBase : public McMove<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMoveBase(McSimulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~ForceBiasMoveBase() = default;

      /**
      * Read body of parameter file block and allocate memory.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before the beginning of each simulation run
      */
      void setup() override;

      /**
      * Attempt and accept or reject a force bias Monte-Carlo move.
      *
      * \return true if accepted, false if rejected
      */
      bool move() override;

      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output() override;

      /**
      * Specify if dc fields need to be saved (returns true).
      */
      bool needsDc() override;

   protected:

      /// Alias for McMove base class.
      using McMoveT = McMove<D,T>;

      // Protected inherited member functions (selected).
      using McMove<D,T>::system;
      using McMove<D,T>::simulator;

   private:

      /// Alias for RField type.
      using RFieldT = typename T::RField;

      /// Local copy of w fields
      DArray< RFieldT > w_;

      /// Copy of initial dc field
      DArray< RFieldT > dc_;

      /// Change in wc
      DArray<RFieldT >  dwc_;

      /// Normal-distributed random field
      RFieldT eta_;

      /// Field used to compute bias for Metropolis criterion
      RFieldT biasField_;

      /// Prefactor of -dc_ in deterministic drift term
      double mobility_;

      /**
      * Compute force bias field.
      */
      virtual
      void computeForceBias(RFieldT& result,
                            RFieldT const & di,
                            RFieldT const & df,
                            RFieldT const & dwc,
                            double mobility) = 0;

   };

   // Public inline methods

   /*
   * Specify if dc fields need to be saved.
   */
   template <int D, class T>
   inline bool ForceBiasMoveBase<D,T>::needsDc()
   {  return true; }

   // Primary template declaration for subclasses
   template <int D, class T> class ForceBiasMove;
}
}
#endif
