#ifndef PSPC_REAL_MOVE_H
#define PSPC_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          //base class
#include <util/global.h>
#include <util/containers/DArray.h>
#include <pspc/field/RField.h>
#include <util/param/ParamComposite.h>


namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * RealMove is a Monte Carlo move in real space
   *
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class RealMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param mcSimulator parent McSimulator
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
      
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
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
      double stepSize_;
      
      /// wField after attempt McMove. local variable wFieldTmp_ used in attemptMove() function
      DArray< RField<D> > wFieldTmp_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
      
   
   };
      

}
}
#endif
