#ifndef FD1D_NR_ITERATOR_H
#define FD1D_NR_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include "Mixture.h"
#include <pscf/LuSolver.h>
#include <util/containers/Array.h>
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

      /**
      * Take one step of Jacobian update.
      */
      void update();

      /**
      * Test criteria for adequate convergence.
      *
      * \return true if converged, false if not converged.
      */
      double residualNorm(Array<double> const & residual) const;

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
