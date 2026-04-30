#ifndef RP_SHIFT_MOVE_H
#define RP_SHIFT_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                    // base class
#include <util/containers/DArray.h>    // member
#include <pscf/math/IntVec.h>          // template with defaults

// Forward declaration
namespace Util {
   template <typename Data> class Array;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * ShiftMove rigidly translates the field.
   * 
   * An attempted ShiftMove rigidly translates all w fields by a random 
   * rigid translation, shifting each coordinate by an integer number of 
   * grid points in each direction.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also both named ShiftMove,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class ShiftMove : public T::McMove
   {

   public:

      // Protected constructor and destructor (see below).

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before the beginning of each simulation run
      */
      void setup() override;

   protected:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMove(typename T::McSimulator& simulator);

      /**
      * Destructor.
      */
      ~ShiftMove() = default;

      /**
      *  Attempt move that translates all w fields.
      */
      void attemptMove() override;

      /**
      * Compute and store shifted w fields.
      *
      * On return, shifted values of fields obtained from system().w()
      * for all monomer types are stored in the w_ member array.
      *
      * \param shift  vector of integer shift values (# of grid points)
      */ 
      virtual void shiftFields(IntVec<D> const & shift) = 0;

      /**
      * Compute a shifted version of a field.
      *
      * This operation is carried out on the CPU.
      * 
      * \param out  shifted field (output)
      * \param in   original field (input)
      * \param shift  vector of integer shift values
      * \param dimensions  dimensions of computational mesh
      */ 
      void shiftField(Array<double>& out, 
                      Array<double> const & in,
                      IntVec<D> shift, 
                      IntVec<D> dimensions) const;

      /// Shifted field configurations.
      mutable DArray< typename T::RField > w_;

      // Alias for McMove base class.
      using McMoveT = typename T::McMove;

      // Inherited protected member functions (selected).
      using McMoveT::system;
      using McMoveT::simulator;

   private:

      /// Maximum absolute value of shift components.
      int maxShift_;

      /// Has private memory been allocated?
      bool isAllocated_;

   };

}
}
#endif
