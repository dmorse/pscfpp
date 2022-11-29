#ifndef PSPC_MIXTURE_H
#define PSPC_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <pscf/inter/Interaction.h>
#include <pscf/chem/Monomer.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>

#include <iostream>

namespace Pscf { 
   template <int D> class Mesh; 
}
 
namespace Pscf {
namespace Pspc
{

   /**
   * Solver for a mixture of polymers and solvents.
   *
   * A Mixture contains a list of Polymer and Solvent objects. Each such
   * object can solve the single-molecule statistical mechanics problem
   * for an ideal gas of the associated species in a set of specified
   * chemical potential fields, and thereby compute concentrations and
   * and single-molecule partition functions. A Mixture is thus both a
   * a chemistry descriptor and an ideal-gas solver.
   *
   * The single-molecule partition functions and concentrations for a
   * non-interacting mixture of polymer and solvent species are computed
   * by invoking the Mixture::compute function.  The Mixture::compute 
   * function takes an arrays of monomer chemical potential fields 
   * (w fields) as an input argument and an array of monomer concentration 
   * fields (c fields) as an output. The objects that store these fields 
   * are owned by the parent System.
   *
   * A Mixture is associated with a Mesh<D> object, which models a spatial
   * discretization mesh, and a UnitCell<D> object, which describes the
   * the periodic unit cell. The Mixture::setupUnitCell function sets up
   * all parameters that depend on the unit cell, and must be called once
   * once after every time the unit cell is initialized or modified, 
   * before the next call to Mixture::compute.
   *
   * \ref user_param_mixture_page "Parameter File Format" 
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Mixture : public MixtureTmpl< Polymer<D>, Solvent<D> >
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
      *
      * This function reads in a complete description of the structure of
      * all species and the composition of the mixture, as well as the
      * target contour length step size ds.
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
      * This function resets unit cell information in the solvers for 
      * every species in the system. It should be called once after
      * every change in the unit cell.
      *
      * \param unitCell UnitCell<D> object that contains Bravais lattice.
      */
      void setupUnitCell(const UnitCell<D>& unitCell);

      /**
      * Reset statistical segment length for one monomer type.
      * 
      * This function resets the kuhn or statistical segment length value
      * for a monomer type, and updates the associcated value in every 
      * block of that monomer type.
      *
      * \param monomerId  monomer type id
      * \param kuhn  new value for the statistical segment length
      */
      void setKuhn(int monomerId, double kuhn);

      /**
      * Compute partition functions and concentrations.
      *
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration fields
      * for blocks of the same monomer type to compute a total monomer
      * concentration (or volume fraction) for each monomer type.
      * Upon return, values are set for volume fraction and chemical 
      * potential (mu) members of each species, and for the 
      * concentration fields for each Block and Solvent. The total
      * concentration for each monomer type is returned in the
      * cFields output parameter. Monomer "concentrations" are returned 
      * in units of inverse steric volume per monomer in an incompressible
      * mixture, and are thus also volume fractions.
      *
      * The arrays wFields and cFields must each have capacity nMonomer(),
      * and contain fields that are indexed by monomer type index. 
      *
      * The optional parameter phiTot is only relevant to problems such as 
      * thin films in which the material is excluded from part of the unit
      * cell by imposing an inhomogeneous constrain on the sum of mononer 
      * concentrations, (i.e., a "mask"). 
      *
      * \param wFields array of chemical potential fields (input)
      * \param cFields array of monomer concentration fields (output)
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(DArray< RField<D> > const & wFields, 
                   DArray< RField<D> >& cFields, 
                   double phiTot = 1.0);
      
      /**
      * Compute derivatives of free energy w/ respect to cell parameters.
      */
      void computeStress();

      /**
      * Combine cFields for each block/solvent into one DArray, which 
      * is used in System.tpp to print a more detailed r-grid file using
      * the command WRITE_C_BLOCK_RGRID.
      * 
      * \param blockCFields empty but allocated DArray to store fields
      */
      void createBlockCRGrid(DArray< RField<D> >& blockCFields) const;

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * Get the pre-computed derivative with respect to unit cell 
      * parameter number n of the free energy per monomer (i.e., of the 
      * product of the free energy density and the monomer reference 
      * volume). The returned value is precomputed by the computeStress()
      * function.
      *
      * \param n  index of unit cell parameter
      */
      double stress(int n) const;

      /**
      * Is this mixture being treated in canonical ensemble?
      *
      * Returns true iff an closed ensemble is used for every polymer
      * and solve species, by specifying a volume fraction phi rather
      * than a chemical potential mu for every species.
      */
      bool isCanonical();

      // Inherited public member functions with non-dependent names
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nMonomer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nPolymer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nSolvent;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nBlock;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::polymer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::monomer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::solvent;

   protected:

      // Inherited protected member functions with non-dependent names
      using MixtureTmpl< Polymer<D>, Solvent<D> >::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Optimal contour length step size.
      double ds_;

      /// Array to store total stress
      FArray<double, 6> stress_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Pointer to associated UnitCell<D>
      UnitCell<D> const * unitCellPtr_;

      /// Return associated domain by reference.
      Mesh<D> const & mesh() const;

      /// Has stress been computed for current w fields?
      bool hasStress_;

   };

   // Inline member function

   // Stress with respect to unit cell parameter n.
   template <int D>
   inline double Mixture<D>::stress(int n) const
   {
      UTIL_CHECK(hasStress_);  
      return stress_[n]; 
   }

   // Get Mesh<D> by constant reference (private).
   template <int D>
   inline Mesh<D> const & Mixture<D>::mesh() const
   {
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   #ifndef PSPC_MIXTURE_TPP
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
