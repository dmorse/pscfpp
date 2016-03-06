#ifndef FD1D_DOMAIN_H
#define FD1D_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include "GeometryMode.h"                  // member

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
      * Generic field type (base class)
      */
      typedef Array<double> Field;

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

      /// \name Accessors
      //@{

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

      /**
      * Get coordinate system flag (Planar, Cylindrical or Spherical).
      */
      GeometryMode const & geometryMode() const;

      //@}
      /// \name Spatial integrals
      //@{

      /**
      * Compute spatial average of a field.
      *
      * \param f a field that depends on one spatial coordinate
      * \return spatial average of field f
      */
      double spatialAverage(Field const & f) const;
 
      /**
      * Compute inner product of two real fields.
      *
      * \param f first field
      * \param g second field
      * \return spatial average of product of two fields.
      */
      double innerProduct(Field const & f, Field const & g) const;

      //@}

   private:

      /**
      * Lower bound of spatial coordinate.
      */
      double xMin_;

      /**
      * Upper bound of spatial coordinate.
      */
      double xMax_;

      /**
      * Spatial discretization step.
      */
      double dx_;

      /**
      * Number of grid points.
      */
      int nx_;

      /**
      * Coordinate system flag (=Planar, Cylindrical, or Spherical).
      */
      GeometryMode geometryMode_;

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

   inline GeometryMode const & Domain::geometryMode() const
   {  return geometryMode_; }

} // namespace Fd1d
} // namespace Pscf
#endif
