#ifndef RP_SHIFT_MOVE_BASE_H
#define RP_SHIFT_MOVE_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>  // base class
#include <util/containers/DArray.h>    // member
#include <pscf/math/IntVec.h>          // template with defaults

// Forward declaration
namespace Util {
   template <typename Data> class Array;
}
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
   }
   namespace Rp {
      template <int D, class T> class System;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * ShiftMove rigidly translates all w fields.
   * 
   * An attempted ShiftMove rigidly translates all w fields by a random 
   * rigid translation, shifting each coordinate by an integer number of 
   * grid points in each direction.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \see \ref rp_ShiftMove_page "Manual Page".
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class ShiftMoveBase : public McMove<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ShiftMoveBase(McSimulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~ShiftMoveBase() = default;

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
      mutable DArray< RField<D,T> > w_;

      // Alias for McMove base class.
      using McMoveT = McMove<D,T>;

      // Inherited protected member functions (selected).
      using McMove<D,T>::system;
      using McMove<D,T>::simulator;

   private:

      /// Maximum absolute value of shift components.
      int maxShift_;

      /// Has private memory been allocated?
      bool isAllocated_;

   };

   // Declaration of primary ShiftMove template
   template <int D, class T> class ShiftMove;

}
}
#endif
