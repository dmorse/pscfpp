#ifndef RP_REAL_MOVE_H
#define RP_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>  // base class
#include <util/containers/DArray.h>    // member
#include <iostream>

// Forward declaration
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * RealMove generates spatially uncorrelated random field changes.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named RealMove, that
   * are defined in Rpc and Rpg namespaces for use in the pscf_rpc and
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : spatial dimension
   *    - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \see \ref rp_RealMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class RealMove : public McMove<D,T>
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      RealMove(McSimulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~RealMove() = default;

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before the simulation loop.
      */
      void setup() override;

   protected:

      /**
      * Attempt unconstrained move.
      */
      void attemptMove() override;

      // Alias for McMove base class.
      using McMoveT = McMove<D,T>;

      // Inherited protected member functions (selected).
      using McMove<D,T>::system;
      using McMove<D,T>::simulator;
      using McMove<D,T>::vecRandom;

   private:

      using RFieldT = RField<D,T>;

      /// New field values, indexed by monomer type.
      DArray< RFieldT > w_;

      /// Change in one field component.
      RFieldT dwc_;

      /// Standard deviation of field changes.
      double sigma_;

      /// Has memory been allocated?
      bool isAllocated_;

   };

}
}
#endif
