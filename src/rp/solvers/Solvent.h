#ifndef RP_SOLVENT_H
#define RP_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/SolventSpecies.h>   // base class
#include <prdc/field/RField.h>          // member
#include <pscf/backends/TmplDeclare.h>  // template declaration macros

// Forward declaration
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Rp {

   using namespace Prdc;

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ref user_param_solvent_sec "Manual Page"
   *
   * <b> Template parameters </b>:
   *
   *   - D : dimension of space (1, 2, or 3)
   *   - T : Types class (CppTp<D> or CudaTp<D>)
   *
   * \ingroup Rp_Solver_Module
   */
   template <int D, class T>
   class Solvent : public SolventSpecies<double>
   {

   public:

      /// Alias for direct (parent) base class.
      using SolventSpeciesT = SolventSpecies<double>;

      /// Alias for indirect (grandparent) base class.
      using SpeciesT = Species<double>;

      // Public member functions

      /**
      * Constructor.
      */
      Solvent();

      /**
      * Destructor.
      */
      ~Solvent();

      /**
      * Create an association with a Mesh.
      *
      * Must be called before allocate().
      *
      * \param mesh  associated Mesh<D> object - spatial discretization
      */
      void associate(Mesh<D> const & mesh);

      /**
      * Allocate memory for a concentration field.
      *
      * This function may only be called once during initialization, after
      * the associate() function and before the first call to the compute()
      * function.
      */
      void allocate();

      /**
      * Compute concentration field, q, and phi or mu.
      *
      * Computes monomer concentration field cField, partition function
      * q, and either the solvent volume fraction phi or solvent chemical
      * potential mu, depending on ensemble. The function takes the
      * chemical potential field wField for the relevant monomer type
      * as its only input argument.
      *
      * The optional parameter phiTot is only relevant to problems such
      * as thin films in which the material is excluded from part of the
      * unit cell by imposing an inhomogeneous constraint on the sum of
      * monomer concentrations (i.e., a "mask"). The default value
      * phiTot = 1.0 is correct for problems with no mask. 
      *
      * The associate and allocate functions must be called before this
      * function may be called.
      *
      * \param wField  monomer chemical potential field of relevant type.
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(RField<D,T> const & wField, double phiTot = 1.0);

      /**
      * Get the monomer concentration field for this solvent.
      *
      * This function retrieves the concentration field computed by the
      * compute() function.
      */
      RField<D,T> const & cField() const;

   private:

      /// Concentration field for this solvent.
      RField<D,T> cField_;

      /// Pointer to associated mesh.
      Mesh<D> const *  meshPtr_;

   };

   /*
   * Get monomer concentration field for this solvent.
   */
   template <int D, class T> inline
   RField<D,T> const & Solvent<D,T>::cField() const
   {  return cField_;  }

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(Solvent)
}
}
#endif
