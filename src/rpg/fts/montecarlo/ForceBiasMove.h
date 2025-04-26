#ifndef RPG_FORCE_BIAS_MOVE_H
#define RPG_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cuda/RField.h>
#include <util/containers/DArray.h> 

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      ForceBiasMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~ForceBiasMove();

      /**
      * Read required parameters from file.
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output();

      /**
      * Setup before the beginning of each simulation run
      */
      void setup();

      /**
      * Attempt and accept or reject force bias Monte-Carlo move.
      *
      * \return true if accepted, false if rejected
      */
      bool move();

      /**
      * Return real move times contributions.
      */
      void outputTimers(std::ostream& out);
      
      /**
      * Decide whether dc fields need to be saved for move
      */
      bool needsDc();

      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;

   protected:

      using McMove<D>::system;
      using McMove<D>::simulator;
      using McMove<D>::random;
      using McMove<D>::cudaRandom;
      using McMove<D>::incrementNAttempt;
      using McMove<D>::incrementNAccept;
      using McMove<D>::incrementNFail;

      using McMove<D>::attemptMoveTimer_;
      using McMove<D>::compressorTimer_;
      using McMove<D>::componentTimer_;
      using McMove<D>::hamiltonianTimer_;
      using McMove<D>::decisionTimer_;
      using McMove<D>::totalTimer_;

   private:

      /// Local copy of w fields
      DArray< RField<D> > w_;

      /// Copy of initial dc field 
      DArray< RField<D> > dc_;

      /// Change in wc
      DArray<RField<D> >  dwc_;

      /// Normal-distributed random fields
      RField<D> gaussianField_;
      
      RField<D> biasField_;
      
      /// Prefactor of -dc_ in deterministic drift term
      double mobility_;
      
   };
   
   // Public inline methods

   /*
   * Return whether dc fields need to be saved for ForceBiasMove.
   */
   template <int D>
   inline bool ForceBiasMove<D>::needsDc()
   {  return true; }

   #ifndef RPG_FORCE_BIAS_MOVE_TPP
   // Suppress implicit instantiation
   extern template class ForceBiasMove<1>;
   extern template class ForceBiasMove<2>;
   extern template class ForceBiasMove<3>;
   #endif

}
}
#endif
