#ifndef FD1D_MIXTURE_H
#define FD1D_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/MixtureTmpl.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Fd1d
{

   class Grid;

   /**
   * Container for species within a mixture.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Mixture : public MixtureTmpl<Polymer, Solvent>
   {

   public:

      // Public typedefs

      /**
      * Monomer chemical potential field type.
      */
      typedef Propagator::WField WField;

      /**
      * Monomer concentration or volume fraction field type.
      */
      typedef Propagator::CField CField;

      // Public member functions

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Destructor.
      */
      ~Mixture();

      /**
      * Read all parameters and initialize.
      *
      * This function reads in a complete description of
      * the chemical composition and structure of all species,
      * as well as the target contour length step size ds.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create an association with the grid and allocate memory.
      * 
      * The grid parameter must have already been initialized, 
      * e.g., by reading its parameters from a file, so that the
      * grid dimensions are known on entry.
      *
      * \param grid associated Grid object (stores address).
      */
      void setGrid(Grid const & grid);

      /**
      * Compute molecular partition functions and concentrations.
      *
      * The arrays wFields and cFields must each have size nMonomer(),
      * and contain fields that are indexed by monomer type index. 
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration
      * fields for blocks of each type to compute a total monomer
      * concentration (or volume fraction) for each monomer type.
      *
      * \param wFields input array of chemical potential fields.
      * \param cFields output array of monomer concentration fields.
      */
      void compute(DArray<WField> const & wFields, DArray<CField>& cFields);

   private:

      /// Optimal contour length step size.
      double ds_;

      /// Pointer to associated Grid object.
      Grid const * gridPtr_;

      /// Return associated grid by reference.
      Grid const & grid() const;

   };

   // Inline member function

   inline Grid const & Mixture::grid() const
   {   
      UTIL_ASSERT(gridPtr_);
      return *gridPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
