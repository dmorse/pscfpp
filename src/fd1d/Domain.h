#ifndef FD1D_GRID_H
#define FD1D_GRID_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * One-dimensional spatial domain and discretization grid.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Domain : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

      /**
      * Set grid parameters. 
      */
      void setParameters(double xMin, double xMax, int nx);

      /**
      * Read all parameters and initialize.
      */
      void readParameters(std::istream& in);

      /**
      * Get minimum spatial coordinate.
      */
      double xMin() const;

      /**
      * Get maximum spatial coordinate.
      */
      double xMax() const;

      /**
      * Get spatial grid step size.
      */
      double dx() const;

      /**
      * Get number of spatial grid points.
      */
      int nx() const;

   private:

      // Lower bound of spatial coordinate
      double xMin_;

      // Upper bound of spatial coordinate
      double xMax_;

      // Spatial discretization step.
      double dx_;

      // Number of grid points.
      int nx_;

   };

   // Inline member functions

   inline int Domain::nx() const
   {  return nx_; }

   inline double Domain::dx() const
   {  return dx_; }

   inline double Domain::xMin() const
   {  return xMin_; }

   inline double Domain::xMax() const
   {  return xMax_; }

} // namespace Fd1d
} // namespace Pscf
#endif
