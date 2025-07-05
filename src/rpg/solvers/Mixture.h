#ifndef RPG_MIXTURE_H
#define RPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/MixtureTmpl.h>     // base class template
#include "Polymer.h"                      // base class parameter
#include "Solvent.h"                      // base class parameter

#include <util/containers/FSArray.h>      // member (stress_)

// Forward declarations
namespace Util { 
   template <typename T> class DArray;
}
namespace Pscf { 
   template <int D> class Mesh; 
   namespace Prdc {
      template <int D> class UnitCell; 
      namespace Cuda {
         template <int D> class WaveList;
         template <int D> class FFT; 
         template <int D> class RField; 
      }
   }
}
 
namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Solver for a mixture of polymers and solvents.
   *
   * A Mixture contain lists of Polymer and Solvent objects. Each such
   * object can solve the single-molecule statistical mechanics problem
   * for an ideal gas of the associated species in a set of specified 
   * chemical potential fields, and thereby compute concentrations and
   * single-molecule partition functions. A Mixture is thus both a 
   * chemistry descriptor and an ideal-gas solver.
   *
   * The single-molecule partition functions and concentrations for a
   * non-interacting mixture of polymer and solvent species are computed
   * by invoking the Mixture::compute function.  The Mixture::compute
   * function takes an arrays of monomer chemical potential fields
   * (w fields) as an input argument and an array of monomer concentration
   * fields (c fields) as an output.
   *
   * A Mixture is associated with a Mesh<D> object, which models a spatial
   * discretization mesh, and a UnitCell<D> object, which models the
   * periodic unit cell. The Mixture::clearUnitCellData function clears
   * parameters that depend on the unit cell as invalid, and must be called
   * once after every time the unit cell parameters are set or modified,
   * before the next call to Mixture::compute.
   *
   * \ingroup Rpg_Solvers_Module
   */
   template <int D>
   class Mixture : public MixtureTmpl< Polymer<D>, Solvent<D> >
   {

   public:

      // Public typename aliases 

      /// Base class.
      using Base = MixtureTmpl< Polymer<D>, Solvent<D> >;

      /// Solvent object type: SolventT = Solvent<D> (inherited)
      using typename Base::SolventT;

      /// Polymer object type: PolymerT = Polymer<D> (inherited)
      using typename Base::PolymerT;

      /// Block type, for a block in a block polymer (inherited)
      using typename Base::BlockT;

      /// Propagator type, for one direction within a block (inherited).
      using typename Base::PropagatorT;

      /// Field type, for data defined on a real-space grid.
      using FieldT = typename PropagatorT::FieldT;

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
      * This function reads in a complete description of the structure of
      * all species and the composition of the mixture, plus a few other 
      * parameters.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create associations with a mesh, FFT, UnitCell, and WaveList object.
      * 
      * The Mesh<D> object must have already been initialized, e.g., by 
      * reading the dimensions from a file, so that the mesh dimensions 
      * are known on entry. The FFT<D> object must have been set up with
      * mesh dimensions equal to those of the mesh. The UnitCell<D> must
      * have been assigned a non-null lattice system, but does not need
      * to have initialized lattice parameters.
      * 
      * The associate and allocate functions are called within the base
      * class readParameters function. This function must be called before 
      * allocate().
      *
      * \param mesh  Mesh<D> object - spatial discretization mesh
      * \param fft  FFT<D> object - Fourier transforms
      * \param cell  UnitCell<D> object - crystallographic unit cell
      * \param waveList  WaveList<D> object - properties of wavevectors
      */
      void associate(Mesh<D> const & mesh, 
                     FFT<D> const & fft, 
                     UnitCell<D> const & cell, 
                     WaveList<D>& waveList);

      /**
      * Allocate required internal memory for all solvers. 
      * 
      * This function is called within the base class readParameters
      * function, after the associate() function.
      */
      void allocate();

      /**
      * Clear data that depends on lattice parameters in all solvers.
      */
      void clearUnitCellData();

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
      * for blocks of each type to compute a total monomer concentrations
      * (or volume fraction) for each monomer type.  Upon return, values 
      * are set for volume fraction (phi) and chemical potential (mu) 
      * members of each species, and for the concentration fields for each 
      * Block and Solvent. The total concentration for each monomer type 
      * is returned in the cFields function parameter. Monomer 
      * "concentrations" fields returned in units of inverse steric volume 
      * per monomer in an incompressible mixture, and are thus also volume 
      * fractions. 
      *
      * The function parameters wFields and cFields must each have array
      * size nMonomer(), and contain fields that are indexed by monomer 
      * type index. 
      * 
      * The optional parameter phiTot is only relevant to problems such as 
      * thin films in which the material is excluded from part of the unit
      * cell by imposing an inhomogeneous constraint on the sum of monomer 
      * concentrations, (i.e., a "mask"). 
      *
      * This function does not compute SCFT free energies or stress (i.e.,
      * derivatives of free energy with respect to unit cell parameters).
      *
      * \param wFields array of chemical potential fields (input)
      * \param cFields array of monomer concentration fields (output)
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(DArray< RField<D> > const & wFields, 
                   DArray< RField<D> > & cFields, 
                   double phiTot = 1.0);

      /**
      * Compute derivatives of free energy w/ respect to cell parameters.
      * 
      * The optional parameter phiTot is only relevant to problems with a 
      * mask, in which the material is excluded from part of the unit cell
      * by imposing an inhomogeneous constrain on the sum of monomer 
      * concentrations. In such cases, the stress needs to be scaled by 
      * a factor of 1/phiTot.
      * 
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void computeStress(double phiTot = 1.0);

      /**
      * Combine cFields for each block/solvent into one DArray.
      *
      * The array created by this function is used by the command
      * WRITE_C_BLOCK_RGRID to write c-fields for all blocks and 
      * species.
      * 
      * \param blockCFields empty but allocated DArray to store fields
      */
      void createBlockCRGrid(DArray< RField<D> >& blockCFields) const;

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
      * Is the ensemble canonical (i.e, closed for all species)?
      *
      * Return true if and only if the species ensemble is closed for all 
      * polymer and solvent species.
      */
      bool isCanonical() const;
	 
      /**
      * Has the stress been computed?
      */
      bool hasStress() const;

      // Public members from MixtureTmpl with non-dependent names 

      using Base::polymer;
      using Base::polymerSpecies;
      using Base::solvent;
      using Base::solventSpecies;
      using MixtureBase::nMonomer;
      using MixtureBase::monomer;
      using MixtureBase::nPolymer;
      using MixtureBase::nSolvent;
      using MixtureBase::nBlock;
      using MixtureBase::vMonomer;

   protected:

      // Inherited protected member functions with non-dependent names 
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      // Private member data

      /// Derivatives of SCFT free energy w/ respect to cell parameters.
      FArray<double, 6> stress_;
   
      /// Target contour length step size (thread model).
      double ds_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Number of unit cell parameters.
      int nParams_;

      /// Has the stress been computed?
      bool hasStress_;

      /// Use batched FFTs to compute stress? (faster, but doubles memory)
      bool useBatchedFFT_;

      // Private member function
      
      /// Return associated Mesh<D> by const reference.
      Mesh<D> const & mesh() const;

   };

   // Inline member function

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

   /*
   * Has the stress been computed?
   */
   template <int D>
   inline bool Mixture<D>::hasStress() const
   {
      return hasStress_;
   }

   #ifndef RPG_MIXTURE_TPP
   // Suppress implicit instantiation
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
//#include "Mixture.tpp"
#endif
