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
      * Compute molecular partition functions and concentrations.
      */
      void compute();

      /**
      * Get number of spatial grid points.
      */
      int nx() const;

      /**
      * Get spatial grid step size.
      */
      double dx() const;

   private:

      // Lower bound of spatial coordinate
      double xMin_;

      // Upper bound of spatial coordinate
      double xMax_;

      // Spatial discretization step.
      double dx_;

      // Optimal contour length step size.
      double ds_;

      // Number of grid points.
      int nx_;

   };

   // Inline member functions

   inline int Mixture::nx() const
   {  return nx_; }

   inline double Mixture::dx() const
   {  return dx_; }

} // namespace Fd1d
} // namespace Pscf
#endif
