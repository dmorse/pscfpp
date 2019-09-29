#ifndef FD1D_SOLVENT_H
#define FD1D_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include <pscf/solvers/SolventTmpl.h>

namespace Pscf { 
namespace Fd1d
{ 

   /**
   * Solver for a point particle solvent species.
   *
   * \ingroup Fd1d_Solver_Module
   */
   class Solvent : public SolventTmpl<Propagator>
   {
   public:

      /**
      * Monomer concentration field type.
      */
      typedef Propagator::CField CField;

      /**
      * Monomer chemical potential field type.
      */
      typedef Propagator::CField WField;

      /**
      * Constructor.
      */
      Solvent();

      /**
      * Destructor.
      */
      ~Solvent();

      void setPhi(double phi);

      void setMu(double mu);

      /**
      * Compute monomer concentration field and partittion function.
      *
      * Upon return, monomer concentration field, phi and mu are set.
      *
      * \param wField monomer chemical potential field
      */
      void compute(WField const & wField);
      //void compute(const DArray<Block::WField>& wFields);

   };

} // namespace Fd1d
} // namespace Pscf
#endif
