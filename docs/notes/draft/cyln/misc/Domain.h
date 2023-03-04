#ifndef CYLN_DOMAIN_H
#define CYLN_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <cyln/field/Field.h>            // member template

namespace Pscf {
namespace Cyln
{

   using namespace Util;

   /**
   * Cylindrical domain and discretization grid.
   *
   * \ingroup Pscf_Cyln_Module
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

      /// \name Initialize parameters
      //@{

      /**
      * Read all parameters and initialize.
      */
      void readParameters(std::istream& in);

      /**
      * Set grid parameters for a planar domain.
      */
      void setParameters(double radius, double length, int nr, int nz);

      //@}
      /// \name Accessors
      //@{

      /**
      * Get radius of cylindrical domain.
      */
      double radius() const;

      /**
      * Get length of domain parallel to rotation axis.
      */
      double length() const;

      /**
      * Get volume of domain.
      */
      double volume() const;

      /**
      * Get radial step size.
      */
      double dr() const;

      /**
      * Get axial step size.
      */
      double dz() const;

      /**
      * Get number of grid points in radial (r) direction.
      */
      int nr() const;

      /**
      * Get number of grid points in axial (z) direction.
      */
      int nz() const;

      //@}
      /// \name Spatial integrals
      //@{

      /**
      * Compute spatial average of a field.
      *
      * \param f a field that depends on one spatial coordinate
      * \return spatial average of field f
      */
      double spatialAverage(Field<double> const & f) const;
 
      /**
      * Compute inner product of two real fields.
      *
      * \param f first field
      * \param g second field
      * \return spatial average of product of two fields.
      */
      double innerProduct(Field<double> const & f, Field<double> const & g) const;

      //@}

   private:

      /**
      * Upper bound of spatial coordinate.
      */
      double radius_;

      /**
      * Lower bound of spatial coordinate.
      */
      double length_;

      /**
      * Generalized D-dimensional volume of simulation cell.
      */
      double volume_;

      /**
      * Radial discretization step.
      */
      double dr_;

      /**
      * Radial discretization step.
      */
      double dz_;

      /**
      * Number of grid points in radial direction.
      */
      int nr_;

      /**
      * Number of grid points in axial (z) direction.
      */
      int nz_;

      /**
      * Total number of grid points.
      */
      int nGrid_;

      /**
      * Work space vector.
      */
      mutable Field<double> work_;

   };

   // Inline member functions

   inline int Domain::nr() const
   {  return nr_; }

   inline int Domain::nz() const
   {  return nz_; }

   inline double Domain::dr() const
   {  return dr_; }

   inline double Domain::dz() const
   {  return dz_; }

   inline double Domain::length() const
   {  return length_; }

   inline double Domain::radius() const
   {  return radius_; }

   inline double Domain::volume() const
   {  return volume_; }

} // namespace Cyln
} // namespace Pscf
#endif
