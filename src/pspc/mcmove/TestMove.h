#ifndef PSPC_TEST_MOVE_H
#define PSPC_TEST_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          //base class
#include <util/global.h>
#include <util/containers/DArray.h>
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class TestMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object.
      */
      TestMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~TestMove();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      void readParameters(std::istream &in);
      
      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output();
      
      // Inherited public member function
      
      using McMove<D>::setup;
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
      
      /// Move step size is randomly selected from uniform distribution [-A, A]
      double A_;
   
   };
      

}
}
#endif
