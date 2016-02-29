#ifndef FD1D_NR_ITERATOR_H
#define FD1D_NR_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include "Mixture.h"
#include <pscf/LuSolver.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * Newton-Raphson Iterator.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class NrIterator : public Iterator
   {

   public:

      /**
      * Monomer chemical potential field.
      */
      typedef Mixture::WField WField;

      /**
      * Monomer concentration / volume fraction field.
      */
      typedef Mixture::CField CField;

      /**
      * Constructor.
      */
      NrIterator();

      /**
      * Destructor.
      */
      virtual ~NrIterator();

      /**
      * Read all parameters and initialize.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Iterate self-consistent field equations to solution.
      * 
      * \return error code: 0 for success, 1 for failure.
      */
      int solve();

      /**
      * Get error tolerance.
      */
      double epsilon();

      /**
      * Compute the residual vector.
      */
      void computeResidual(Array<WField> const & wFields, 
                           Array<WField> const & cFields, 
                           Array<double>& residual);

      void computeJacobian();

      void update();

      /**
      * Test criteria for adequate convergence.
      *
      * \return true if converged, false if not converged.
      */
      bool isConverged();

   private:

      LuSolver solver_;

      /// Perturbed chemical potential fields
      DArray<WField> wFieldsNew_;

      /// Perturbed monomer concentration fields
      DArray<WField> cFieldsNew_;

      DArray<double> cArray_;
      DArray<double> wArray_;

      /// Residual vector. size = (# monomers)x(# grid points).
      DArray<double> residual_;

      /// Jacobian matrix
      DMatrix<double> jacobian_;

      /// Perturbed residual
      DArray<double> residualNew_;

      /// Change in field
      DArray<double> dOmega_;

      /// Error tolerance.
      double epsilon_;

      /// Have arrays been allocated?
      bool isAllocated_;

      /**
      * Allocate memory if needed. If isAllocated, check array sizes.
      */
      void allocate();

   };

   // Inline function

   inline double NrIterator::epsilon()
   {  return epsilon_; }

} // namespace Fd1d
} // namespace Pscf
#endif
