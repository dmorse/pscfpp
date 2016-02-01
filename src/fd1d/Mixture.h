#ifndef FD1D_MIXTURE_H
#define FD1D_MIXTURE_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/MixtureTmpl.h>

namespace Pscf {
namespace Fd1d
{

   class Grid;

   class Mixture : public MixtureTmpl<Polymer, Solvent>
   {

   public:

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
      */
      void readParameters(std::istream& in);

      /**
      * Set grid and allocate all required memory.
      * 
      * \param grid associated Grid object (stores address).
      */
      void setGrid(Grid const & grid);

      /**
      * Compute molecular partition functions and concentrations.
      */
      void compute();

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
