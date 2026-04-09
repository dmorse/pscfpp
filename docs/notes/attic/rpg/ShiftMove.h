#ifndef RPG_SHIFT_MOVE_H
#define RPG_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cpu/RField.h>                 // member
#include <util/containers/DArray.h>          // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * ShiftMove shifts field.
   *
   * \see \ref rp_ShiftMove_page "Manual Page". 
   *
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class ShiftMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      ShiftMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~ShiftMove();

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
      * Return field shift move times contributions.
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
      using McMove<D>::simulator;
      using McMove<D>::random;

      /**
      *  Attempt shift field move.
      *
      *  This function should shift the system w fields in r-grid
      *  format, as returned by system().w().rgrid()
      *
      */
      void attemptMove();

   private:
      
      // Initial field values
      DArray< RField<D> > w0_;

      // New field values after shift
      DArray< RField<D> > w_;
      
      // The shift range
      int maxShift_;
      
      // Has the variable been allocated?
      bool isAllocated_;
   
   };
      
   // Explicit instantiation declarations
   extern template class ShiftMove<1>;
   extern template class ShiftMove<2>;
   extern template class ShiftMove<3>;

}
}
#endif
