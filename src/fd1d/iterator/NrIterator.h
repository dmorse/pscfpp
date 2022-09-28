#ifndef FD1D_NR_ITERATOR_H
#define FD1D_NR_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <fd1d/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/Array.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * Newton-Raphson Iterator for SCF equations.
   *
   * \ingroup Fd1d_Iterator_Module
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
      * Default constructor.
      */
      NrIterator();

      /**
      * Constructor.
      *
      * \param system parent System object.
      */
      NrIterator(System& system);

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
      * \param isContinuation True if part of sweep, and not first step.
      * \return error code: 0 for success, 1 for failure.
      */
      int solve(bool isContinuation = false);

      /**
      * Get error tolerance.
      */
      double epsilon();

      /**
      * Compute the residual vector.
      *
      * \param wFields monomer chemical potential fields (input)
      * \param cFields monomer concentration fields (input)
      * \param residual vector of residuals (errors) (output)
      */
      void computeResidual(Array<WField> const & wFields, 
                           Array<WField> const & cFields, 
                           Array<double>& residual);

      /**
      * Compute and return norm of a residual vector.
      *
      * \param residual vector of residuals (errors) (input)
      */
      double residualNorm(Array<double> const & residual) const;

      /**
      * Compute the Jacobian matrix (stored in class member).
      */
      void computeJacobian();

   private:

      /// Solver for linear system Ax = b.
      LuSolver solver_;

      /// Perturbed chemical potential fields (work space).
      DArray<WField> wFieldsNew_;

      /// Perturbed monomer concentration fields (work space).
      DArray<WField> cFieldsNew_;

      /// Concentrations at one point (work space).
      DArray<double> cArray_;

      /// Chemical potentials at one point (work space).
      DArray<double> wArray_;

      /// Residual vector. size = nr = (# monomers)x(# grid points).
      DArray<double> residual_;

      /// Jacobian matrix. Dimensions nr x nr.
      DMatrix<double> jacobian_;

      /// Perturbed residual. size = nr.
      DArray<double> residualNew_;

      /// Change in field
      DArray<double> dOmega_;

      /// Error tolerance.
      double epsilon_;

      /// Have arrays been allocated?
      bool isAllocated_;

      /// Has the Jacobian been calculated at current state?
      bool newJacobian_;

      /// Does the Jacobian need to be re-calculated ?
      bool needsJacobian_;

      /// Is the ensemble canonical for all species ?
      bool isCanonical_;

      /**
      * Allocate memory if needed. If isAllocated, check array sizes.
      */
      void allocate();

      /**
      * Increment the chemical potential fields
      *
      * \param wOld array of old chemical potential fields
      * \param dW array of increments, indexed as in residual columns
      * \param wNew array of new chemical potential fields
      */
      void incrementWFields(Array<WField> const & wOld,
                            Array<double> const & dW,
                            Array<WField>& wNew);

   };

   // Inline function

   inline double NrIterator::epsilon()
   {  return epsilon_; }

} // namespace Fd1d
} // namespace Pscf
#endif
