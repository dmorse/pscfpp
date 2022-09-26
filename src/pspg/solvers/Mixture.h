#ifndef PSPG_MIXTURE_H
#define PSPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/DArray.h>

namespace Pscf { 
   template <int D> class Mesh; 
}
 
namespace Pscf {
namespace Pspg
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
   * \ingroup Pspg_Solvers_Module
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
      * \param unitCell crystallographic unit cell
      * \param waveList container for wavevector data
      */
      void setupUnitCell(UnitCell<D> const & unitCell, 
                         WaveList<D> const & waveList);

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
      compute(DArray< RDField<D> > const & wFields, 
              DArray< RDField<D> > & cFields);

      /**
      * Get monomer reference volume.
      *
      * \param waveList container for wavevector data
      */
      void computeStress(WaveList<D> const & waveList);

      /**
      * Combine cFields for each block/solvent into one DArray.
      *
      * The array created by this function is used by the command
      * WRITE_C_BLOCK_RGRID to write c-fields for all blocks and 
      * species.
      * 
      * \param blockCFields empty but allocated DArray to store fields
      */
      void createBlockCRGrid(DArray< RDField<D> >& blockCFields) const;

      /**
      * Get derivative of free energy w/ respect to cell parameter.
      *
      * Get precomputed value of derivative of free energy per monomer
      * with respect to unit cell parameter number n.
      *
      * \param parameterId  unit cell parameter index
      */
      double stress(int parameterId) const;

      /**
      * Get monomer reference volume.
      */
      double vMonomer() const;

      /**
      * Is the ensemble canonical (i.e, closed for all species)?
      *
      * Return true if and only if the ensemble is closed for all polymer 
      * and solvent species.
      */
      bool isCanonical();
	 
      // Public members from MixtureTmpl with non-dependent names 
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nMonomer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nPolymer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nSolvent;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::nBlock;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::polymer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::monomer;
      using MixtureTmpl< Polymer<D>, Solvent<D> >::solvent;

   protected:

      // Public members from MixtureTmpl with non-dependent names 
      using MixtureTmpl< Polymer<D>, Solvent<D> >::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Derivatives of free energy w/ respect to cell parameters.
      FArray<double, 6> stress_;
    
      /// Monomer reference volume (set to 1.0 by default).
      double vMonomer_;

      /// Optimal contour length step size.
      double ds_;

      /// Number of unit cell parameters.
      int nUnitCellParams_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Has the stress been computed?
      bool hasStress_;

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

   /*
   * Get derivative of free energy w/ respect to cell parameter.
   */
   template <int D>
   inline double Mixture<D>::stress(int parameterId) const
   {
      UTIL_CHECK(hasStress_);  
      return stress_[parameterId]; 
   }

   #ifndef PSPG_MIXTURE_TPP
   // Suppress implicit instantiation
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
//#include "Mixture.tpp"
#endif
