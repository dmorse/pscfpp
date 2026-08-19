#ifndef R1D_HOMOGENEOUS_COMPARISON_H
#define R1D_HOMOGENEOUS_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <r1d/system/SystemAccess.h>             // base class

namespace Pscf {
   class FhMixture;
   class FhInteraction;
}

namespace Pscf {
namespace R1d {

   /**
   * Command to compute properties of homogeneous reference system.
   *
   * \ingroup R1d_System_Module
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
      * Compute computed properties of a homogeneous reference system.
      */
      void clear();

      /**
      * Output comparison to a homogeneous reference system.
      *
      * \param mode mode index
      * \param out output stream 
      */
      void output(int mode, std::ostream& out) const;

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

      /**
      * Pointer to a FhMixture.
      */
      FhMixture* fhMixturePtr_;

      /**
      * Pointer to a private FhInteraction.
      */
      FhInteraction* fhInteractionPtr_;

      /**
      * Get the FhMixture by reference.
      */
      FhMixture& fhMixture() { return *fhMixturePtr_; }

      /**
      * Get the FhMixture by const reference.
      */
      FhMixture const & fhMixture() const 
      { return *fhMixturePtr_; }

      /**
      * Get the FhInteraction by reference.
      */
      FhInteraction& fhInteraction() { return *fhInteractionPtr_; }

      /**
      * Get the FhInteraction by const reference.
      */
      FhInteraction const & fhInteraction() const
      { return *fhInteractionPtr_; } 

      bool hasData_;

   };

} // namespace R1d
} // namespace Pscf
#endif
