#ifndef FD1D_GRID_H
#define FD1D_GRID_H

/*
* PFTS - Polymer Field Theory Simulator
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
   * Spatial discretization grid for a one-dimensional problem.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Grid : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Grid();

      /**
      * Destructor.
      */
      ~Grid();

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

   inline int Grid::nx() const
   {  return nx_; }

   inline double Grid::dx() const
   {  return dx_; }

   inline double Grid::xMin() const
   {  return xMin_; }

   inline double Grid::xMax() const
   {  return xMax_; }

} // namespace Fd1d
} // namespace Pscf
#endif
