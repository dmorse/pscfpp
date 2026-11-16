#ifndef R1_POLYMER_H
#define R1_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PolymerTmpl.h>   // base class template
#include <util/containers/DArray.h>      // member template

// Forward declarations
namespace Pscf {
   namespace R1 {
      class Block;
      class Propagator;
   }
}

namespace Pscf {
namespace R1 {

   using namespace Util;

   /**
   * Descriptor and solver for a block polymer species.
   *
   * A Polymer contains an array of Block objects, each of which contains
   * two propagator objects. After Polymer::compute is called, each block
   * stores the block concentration (volume fraction) fields computed for 
   * that block, and each Propagator contains the solution of the modified
   * diffusion equation for that propagator. Functions that return values
   * for molecular volume fraction phi, molecular chemical potential mu, 
   * and molecular partition function q are indirectly inherited from 
   * Pscf::PolymerSpecies.
   * 
   * \ref user_param_polymer_sec "Parameter File Format"
   * \ingroup R1_Solver_Module
   */
   class Polymer : public PolymerTmpl<Block, Propagator, double>
   {

   public:

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Compute solution to modified diffusion equation and concentrations.
      *
      * Upon return, propagator solutions and block concentration fields
      * have been computed for all propagators and blocks, along with the
      * molecular partition function q and phi or mu.
      *
      * Each block solves the MDE in the field w of its monomer type
      * minus a constant shift (see Block::setupSolver), so the stored
      * propagators are those of the shifted fields. Their partition
      * function is qShifted(); ln q, mu, phi and the block
      * concentrations do not depend on the shifts. The shift is
      *
      *    min(w) - ln < exp(-N (w - min(w))) > / N,
      *
      * where < > is the spatial (volume) average and N is length(): the
      * free energy per monomer of a molecule of length N made only of
      * monomers of that type, held at one point. It lies between the
      * minimum and the average of w, so that the propagators stay
      * bounded and qShifted() is at most about one. It follows a
      * constant added to w exactly, and it is close to the field where
      * the molecules are, which keeps the Crank-Nicolson contour error
      * small. If qShifted() is nevertheless below about 1e-294 (a long
      * chain in a field that varies strongly over distances that it
      * cannot follow), the shifts are raised so that qShifted() is one
      * (see PolymerTmpl::raiseShifts). A raised shift that exceeds the
      * field by more than 2/ds at some grid point makes the
      * Crank-Nicolson factor negative there, which requires a field
      * that varies by more than 2/ds.
      *
      * \param wFields  array of chemica potential fields.
      */
      void compute(DArray< DArray<double> > const & wFields);

   private:

      /// Work array for the field shifts (allocated by compute).
      DArray<double> wWork_;

      // Restrict access
      using PolymerTmpl<Block, Propagator, double>::solve;

   };

} // namespace R1
} // namespace Pscf

// Explicit instantiation declaration
namespace Pscf {
   extern template class PolymerTmpl<R1::Block, R1::Propagator, double>;
}
#endif
