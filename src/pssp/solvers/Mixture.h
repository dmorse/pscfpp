#ifndef PSSP_MIXTURE_H
#define PSSP_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>

namespace Pscf { 
   template <int D> class Mesh; 
   namespace Pssp{
      template <int D> class Basis;
   }
}
 
namespace Pscf {
namespace Pssp
{


   /**
   * Solver for a mixture of polymers and solvents.
   *
   * A Mixture contains a list of Polymer and Solvent objects. Each
   * such object can solve the single-molecule statistical mechanics 
   * problem for an ideal gas of the associated species in a set of
   * specified chemical potential fields, and thereby compute 
   * concentrations and single-molecule partition functions. A
   * Mixture is thus both a chemistry descriptor and an ideal-gas 
   * solver.
   *
   * A Mixture is associated with a Mesh<D> object, which models a
   * spatial discretization mesh. 
   *
   * \ingroup Pssp_Solvers_Module
   */
   template <int D>
   class Mixture : public MixtureTmpl< Polymer<D>, Solvent<D> >
   {

   public:

      // Public typedefs

      /**
      * Monomer chemical potential field type.
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Monomer concentration or volume fraction field type.
      */
      typedef typename Propagator<D>::CField CField;

      // Public member functions

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
      *
      * This function reads in a complete description of
      * the chemical composition and structure of all species,
      * as well as the target contour length step size ds.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create an association with the mesh and allocate memory.
      * 
      * The Mesh<D> object must have already been initialized, 
      * e.g., by reading its parameters from a file, so that the
      * mesh dimensions are known on entry.
      *
      * \param mesh associated Mesh<D> object (stores address).
      */
      void setMesh(Mesh<D> const & mesh);

      /**
      * Set unit cell parameters used in solver.
      * 
      * \param unitCell UnitCell<D> object that contains Bravais lattice.
      */
      void setupUnitCell(const UnitCell<D>& unitCell);

      /**
      * Compute concentrations.
      *
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration
      * fields for blocks of each type to compute a total monomer
      * concentration (or volume fraction) for each monomer type.
      * Upon return, values are set for volume fraction and chemical 
      * potential (mu) members of each species, and for the 
      * concentration fields for each Block and Solvent. The total
      * concentration for each monomer type is returned in the
      * cFields output parameter. 
      *
      * The arrays wFields and cFields must each have size nMonomer(),
      * and contain fields that are indexed by monomer type index. 
      *
      * \param wFields array of chemical potential fields (input)
      * \param cFields array of monomer concentration fields (output)
      */
      void 
      compute(DArray<WField> const & wFields, DArray<CField>& cFields);
      
      // Array to store total stress
      FArray<double, 6> TStress;

      // Function to calculate total stress of the unit cell
      void
      computeTStress(Basis<D>& basis);

      /**
      * Get monomer reference volume.
      */
      double vMonomer() const;

      using MixtureTmpl< Polymer<D>, Solvent<D> >::nMonomer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nPolymer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nSolvent;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::polymer;

   protected:

      using MixtureTmpl< Polymer<D>, Solvent<D> >::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Monomer reference volume (set to 1.0 by default).
      double vMonomer_;

      /// Optimal contour length step size.
      double ds_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Return associated domain by reference.
      Mesh<D> const & mesh() const;

   };

   // Inline member function

   /*
   * Get monomer reference volume (public).
   */
   template <int D>
   inline double Mixture<D>::vMonomer() const
   {  return vMonomer_; }

   /*
   * Get Mesh<D> by constant reference (private).
   */
   template <int D>
   inline Mesh<D> const & Mixture<D>::mesh() const
   {   
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   #ifndef PSSP_MIXTURE_TPP
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;
   #endif

} // namespace Pssp
} // namespace Pscf
// #include "Mixture.tpp"
#endif
