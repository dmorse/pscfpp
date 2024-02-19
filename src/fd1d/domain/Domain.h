#ifndef R1D_DOMAIN_H
#define R1D_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include "GeometryMode.h"                  // member
#include <util/containers/DArray.h>        // function parameter

namespace Pscf {
namespace R1d {

   using namespace Util;

   /**
   * One-dimensional spatial domain and discretization grid.
   *
   * \ref fd1d_Domain_page "Parameter File Format"
   * \ingroup R1d_Domain_Module
   */
   class Domain : public ParamComposite
   {

   public:

      /**
      * Generic field type (base class)
      */
      typedef DArray<double> Field;

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

      /// \name Initialize parameters
      ///@{

      /**
      * Read all parameters and initialize.
      */
      void readParameters(std::istream& in);

      /**
      * Set grid parameters for a planar domain.
      *
      * \param xMin  minimum normal coordinate value
      * \param xMax  maximum normal coordinate value
      * \param nx  number of grid points, including endpoints
      */
      void setPlanarParameters(double xMin, double xMax, int nx);

      /**
      * Set grid parameters for a cylindrical or spherical shell.
      *
      * \param mode  enumeration (Llamellar, Cylindrical or Spherical)
      * \param xMin  minimum radius
      * \param xMax  maximum radius
      * \param nx  number of grid points, including endpoints
      */
      void setShellParameters(GeometryMode mode, 
                              double xMin, double xMax, 
                              int nx);

      /**
      * Set grid parameters for a cylinder.
      * 
      * \param xMax  maximum radius
      * \param nx  number of grid points, including endpoints
      */
      void setCylinderParameters(double xMax, int nx);

      /**
      * Set grid parameters for a sphere.
      *
      * \param xMax  maximum radius
      * \param nx  number of grid points, including endpoints
      */
      void setSphereParameters(double xMax, int nx);

      ///@}
      /// \name Accessors
      ///@{

      /**
      * Get minimum spatial coordinate.
      */
      double xMin() const;

      /**
      * Get maximum spatial coordinate.
      */
      double xMax() const;

      /**
      * Get generalized volume of domain.
      *
      * Returns volume of spherical domain, area of cylindrical 
      * domain, or a length of a planar domain.
      */
      double volume() const;

      /**
      * Get spatial grid step size.
      */
      double dx() const;

      /**
      * Get number of spatial grid points, including both endpoints.
      */
      int nx() const;

      /**
      * Get coordinate system flag (Planar, Cylindrical or Spherical).
      */
      GeometryMode const & mode() const;

      /**
      * Is this a cylindrical or spherical shell?
      *
      * This value is relevant only if the geometry mode is spherical or
      * cylindrical. If so, isShell is set true if the optional parameter
      * xMin is present and assigned a positive value in the parameter
      * file. If geometryMode is planar or xMin is absent, then isShell
      * is false.
      */
      bool isShell() const;

      ///@}
      /// \name Spatial integrals
      ///@{

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

      ///@}

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
      * Generalized D-dimensional volume of simulation cell.
      */
      double volume_;

      /**
      * Number of grid points.
      */
      int nx_;

      /**
      * Coordinate system flag (=Planar, Cylindrical, or Spherical).
      */
      GeometryMode mode_;

      /**
      * Is this a cylindrical or spherical shell?
      */
      bool isShell_;

      /**
      * Work space vector.
      */
      mutable DArray<double> work_;

      /**
      * Compute generalized volume, called by each set function.
      */
      void computeVolume();

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

   inline double Domain::volume() const
   {  return volume_; }

   inline GeometryMode const & Domain::mode() const
   {  return mode_; }

   inline bool Domain::isShell() const
   {  return isShell_; }

} // namespace R1d
} // namespace Pscf
#endif
