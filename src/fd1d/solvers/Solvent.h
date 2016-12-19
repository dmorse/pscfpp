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
   * \ingroup Pscf_Fd1d_Module
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

      /**
      * Compute monomer concentration field and partittion function.
      *
      * Upon return, monomer concentration field, phi and mu are set.
      *
      * \param wField monomer chemical potential field
      */
      void compute(WField const & wField);

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const
      {  return concentration_;  }
   
   private:
   
      CField concentration_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
