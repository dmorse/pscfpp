#ifndef PSPG_REAL_MOVE_H
#define PSPG_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          //base class
#include <util/global.h>
#include <util/param/ParamComposite.h>
#include <pspg/field/RField.h>
#include <pspg/field/Field.h>  
#include <util/containers/DArray.h>
#include <curand.h>


namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /**
   * RealMove is a Monte Carlo move in real space
   *
   * \ingroup Pspg_Simulate_McMove_Module
   */
   template <int D>
   class RealMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      RealMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~RealMove();

      /**
      * Read required parameters from file.
      *
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
      * Return real move times contributions.
      */
      void outputTimers(std::ostream& out);
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;
      

   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;
      /**
      *  Attempt unconstrained move.
      *
      *  This function should modify the system w fields in r-grid
      *  format, as returned by system().w().rgrid(), in order apply
      *  an unconstrained attempted move. The compressor will then be
      *  applied in order to restore the density constraint.
      *
      */
      void attemptMove();


   private:
      
      /// Move step size is randomly selected from uniform distribution [-stepSize_, stepSize_]
      float stepSize_;
      
      /// Random fields between [-stepSize_, stepSize_]
      RField<D> randomField_;
      
      /// wField after attempt McMove. local variable wFieldTmp_ used in attemptMove() function
      DArray< RField<D> > wFieldTmp_;
      
      /// Has the variable been allocated?
      bool isAllocated_;
      
      /// GPU random number generator
      curandGenerator_t gen_;
      
   
   };
      

}
}
#endif
