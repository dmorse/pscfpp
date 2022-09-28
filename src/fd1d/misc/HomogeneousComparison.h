#ifndef FD1D_HOMOGENEOUS_COMPARISON_H
#define FD1D_HOMOGENEOUS_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fd1d/SystemAccess.h>             // base class

namespace Pscf {
namespace Fd1d {

   /**
   * Command to compute properties of homogeneous reference system.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class HomogeneousComparison : public SystemAccess
   {

   public:

      /**
      * Default constructor.
      */
      HomogeneousComparison();

      /**
      * Constructor.
      */
      HomogeneousComparison(System& system);

      /**
      * Destructor.
      */
      ~HomogeneousComparison();

      /**
      * Compute properties of a homogeneous reference system.
      *
      * This function should be called after iterator().solve()
      * to compute properties of a homogeneous reference system
      * to which the properties of the system of interest can 
      * be compared. The value of the mode parameter controls
      * the choice of homogeneous reference system used for this
      * comparison.
      *
      * Mode parameter values:
      *
      *    - mode = 0   : homogeneous system with same phi's
      *    - mode = 1,2 : homogeneous system with same mu's
      *
      * The difference between mode indices 1 and 2 is the 
      * initial guess used in the iterative computation of
      * the composition of the homogeneous reference system:
      *
      *    - mode = 1  : composition at last grid point (nx -1)
      *    - mode = 2  : composition at first grid point (0)
      * 
      * Mode indices 1 and 2 are intended to be used for 
      * calculation of excess properties in, e.g., computation
      * of properties of a micelle or an interface.
      *
      * \param mode mode index
      */
      void compute(int mode);

      /**
      * Output comparison to a homogeneous reference system.
      *
      * \param mode mode index
      * \param out output stream 
      */
      void output(int mode, std::ostream& out);

   private:

      /**
      * Work array (size = # of grid points).
      */
      DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> c_;

      /**
      * Work array (size = # of molecular species).
      */
      DArray<double> p_;

      /**
      * Work array (size = # of molecular species).
      */
      DArray<double> m_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
