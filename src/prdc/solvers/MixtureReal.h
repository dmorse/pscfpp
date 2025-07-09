#ifndef PRDC_MIXTURE_REAL_H
#define PRDC_MIXTURE_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/MixtureTmpl.h>     // base class template
#include <util/containers/FArray.h>       // member template (stress_)
#include <iostream>

// Forward declarations
namespace Util {
   template <typename T> class DArray;
}
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Solver and descriptor for a mixture of polymers and solvents.
   *
   * A MixtureReal contains lists of Polymer (PT) and Solvent (ST) 
   * objects. Each such object can solve statistical mechanics of a single 
   * molecule of the associated species in a set of specified chemical 
   * potential fields, and thereby compute concentrations and molecular 
   * partition functions for all species in non-interacting reference 
   * system. 
   *
   * The compute() member function computes single-molecule partition 
   * functions and monomer concentrations all species.  This function 
   * takes an array of chemical potential fields (w fields) acting on
   * different monomer types an  input and yields an array of total
   * monomer concentration fields (c fields) as an output.
   *
   * \ref user_param_mixture_page "Manual Page"
   */
   template <int D, class PT, class ST>
   class MixtureReal : public MixtureTmpl<PT, ST>
   {

   public:

      // Public type name aliases

      /// MixtureTmplT class.
      using MixtureTmplT = MixtureTmpl<PT,ST>;

      /// Solvent object type: SolventT = ST (inherited).
      using typename MixtureTmplT::SolventT;

      /// Polymer object type: PolymerT = PT (inherited).
      using typename MixtureTmplT::PolymerT;

      /// Block type, for a block in a block polymer (inherited).
      using typename MixtureTmplT::BlockT;

      /// Propagator type, for one direction within a block (inherited).
      using typename MixtureTmplT::PropagatorT;

      /// Field type, for data defined on a real-space grid.
      using FieldT = typename PropagatorT::FieldT;

      /// WaveList type.
      using FFTT = typename BlockT::FFTT;

      /// WaveList type.
      using WaveListT = typename BlockT::WaveListT;

      /// FieldIo type.
      using FieldIoT = typename BlockT::FieldIoT;

      // Public member functions

      /// \name Construction, Initialization and Destruction
      ///@{

      /**
      * Constructor.
      */
      MixtureReal();

      /**
      * Destructor.
      */
      ~MixtureReal();

      /**
      * Read all parameters and initialize.
      *
      * This function reads in a complete description of the structure of
      * all species and the composition of the mixture, plus a few other
      * parameters.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Create associations with Mesh, FFT, UnitCell, and WaveList objects.
      *
      * The Mesh<D> object must have already been initialized, e.g., by
      * reading the dimensions from a file, so that the mesh dimensions
      * are known on entry. The FFTT object must have been set up with
      * mesh dimensions equal to those of the mesh. The UnitCell<D> must
      * have been assigned a non-null lattice system, but does not need
      * to have initialized lattice parameters.
      *
      * This function is called within the readParameters function of the
      * parent System, after calls to the readParameters member functions 
      * of the Mixture and Domain. The Mesh, FFT, and UnitCell lattice are
      * initialized by the Domain readParameters function prior to entry.
      *
      * \param mesh  associated Mesh<D> object
      * \param fft  associated FFTT object (Fast Fourier Transform type)
      * \param cell  associated UnitCell<D> object
      * \param waveList  associated WaveListT object
      */
      void associate(Mesh<D> const & mesh,
                     FFTT const & fft,
                     UnitCell<D> const & cell,
                     WaveListT& waveList);

      /**
      * Create an association with a FieldIoT object.
      *
      * The associated FieldIoT is only used by member functions that write 
      * concentration or propagator fields to file.
      *
      * \param fieldIo  associated FieldIoT object
      */
      void setFieldIo(FieldIoT const & fieldIo);

      /**
      * Allocate required internal memory for all solvers.
      *
      * This function is called within the readParameters of the parent
      * System, after the associate() function.
      */
      void allocate();

      ///@}
      /// \name Primary Computations
      ///@{

      /**
      * Compute partition functions and concentrations.
      *
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration fields
      * for blocks of each type to compute a total monomer concentration
      * (or volume fraction) for each monomer type.  Upon return, values
      * are set for volume fraction (phi) and chemical potential (mu)
      * members of each species, and for the concentration fields for each
      * Block and Solvent. The total concentration for each monomer type is
      * returned in the cFields function parameter. Monomer concentration
      * fields are normalized by the inverse steric volume per monomer in
      * an incompressible mixture, and are thus also volume fractions.
      *
      * The array function parameters wFields and cFields must each have
      * capacity nMonomer(), and contain fields that are indexed by monomer
      * type index.
      *
      * The optional parameter phiTot is only relevant to problems such as
      * thin films in which the material is excluded from part of the unit
      * cell by imposing an inhomogeneous constraint on the sum of monomer
      * concentrations, (i.e., a "mask"). In such cases, the volume
      * fraction phi for each species is defined as a fraction of the
      * volume that is occupied by material, rather than as a fraction
      * of the entire unit cell volume.
      *
      * This function does not compute SCFT free energies or stress (i.e.,
      * derivatives of free energy with respect to unit cell parameters).
      *
      * \param wFields  array of chemical potential fields (input)
      * \param cFields  array of monomer concentration fields (output)
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(DArray<FieldT> const & wFields,
                   DArray<FieldT>& cFields,
                   double phiTot = 1.0);

      /**
      * Set the isSymmetric flag true or false.
      *
      * The isSymmetric variable affects whether a space group is written
      * in field file headers by functions that write concentration or
      * propagator fields to file.  This variable should be set true after
      * calling compute if the w fields passed to the compute function 
      * were known to be symmetric under the space group. 
      */
      void setIsSymmetric(bool isSymmetric);

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
      * Has the stress been computed since the last MDE solution?
      */
      bool hasStress() const;

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * Get the pre-computed derivative of the free energy per monomer with
      * respect to a single unit cell parameter.
      *
      * \param parameterId  index of unit cell parameter
      */
      double stress(int parameterId) const;

      ///@}
      /// \name Parameter Modification
      ///@{

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
      * Clear all data that depends on the unit cell parameters.
      *
      * This function marks all private data that depends on the values of
      * the unit cell parameters as invalid, so that it can be recomputed
      * before it is next needed. This function should be called after any 
      * change in unit cell parameters.
      */
      void clearUnitCellData();

      ///@}
      /// \name Concentration Field Output
      ///@{

      /**
      * Get c-fields for all blocks and solvents as array of r-grid fields.
      *
      * On return, each element of the blockCFields array contains the
      * monomer concentration field for a single block of a polymer species
      * or a single solvent species. These are indexed with polymer blocks
      * first, followed by solvent species. Polymer blocks are listed with
      * blocks of each polymer placed consecutively in order of block
      * index, with polymers ordered by polymer index. Fields associated
      * with solvents are listed after all polymer blocks, ordered by
      * solvent species index.
      *
      * This function will allocate the blockCFields array and the
      * FieldT arrays it contains as needed. This array thus does not
      * need to be allocated on entry. If the array or the fields objects
      * it contains are allocated on entry, their capacities must be
      * correct or an error will be thrown.
      *
      * \param blockCFields  DArray of FieldT field objects (output)
      */
      void createBlockCRGrid(DArray<FieldT>& blockCFields) const;

      /**
      * Write c fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with polymer
      * blocks appear first, followed by columns associated with solvent
      * species. Polymer blocks are listed with blocks blocks of each 
      * species ordered consecutively in order of block index, with groups
      * of columns associated with polymers ordered by polymer index. 
      * Solvents are ordered by solvent species index.
      *
      * \param filename  name of output file
      */
      void writeBlockCRGrid(std::string const & filename) const;

      ///@}
      /// \name Propagator Output
      ///@{
      
      /**
      * Write one slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write the complete propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name. The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

      ///@}

      // Inherited public member functions
      using MixtureTmplT::polymer;
      using MixtureTmplT::polymerSpecies;
      using MixtureTmplT::solvent;
      using MixtureTmplT::solventSpecies;
      using MixtureBase::nMonomer;
      using MixtureBase::monomer;
      using MixtureBase::nPolymer;
      using MixtureBase::nSolvent;
      using MixtureBase::nBlock;
      using MixtureBase::vMonomer;
      using MixtureBase::isCanonical;

   protected:

      /// Return associated Mesh<D> by const reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

      /// Return associated UnitCell<D> by const reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Return associated FieldIoT by const reference.
      FieldIoT const & fieldIo() const
      {  return *fieldIoPtr_; }

      /// Return target value for the contour step size ds.
      double ds() const
      {  return ds_; }

      // Inherited protected member functions
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      // Private member data

      /// Derivatives of SCFT free energy w/ respect to cell parameters.
      FArray<double, 6> stress_;

      /// Target contour length step size (only used in thread model).
      double ds_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      // Pointer to associated UnitCell<D> object 
      UnitCell<D> const * unitCellPtr_;

      // Pointer to associated FieldIoT object 
      FieldIoT const * fieldIoPtr_;

      /// Number of unit cell parameters.
      int nParam_;

      /// Has stress been computed for the current w fields?
      bool hasStress_;

      // Set true iff the w fields used in the MDE are symmetric
      bool isSymmetric_;

      // Private member functions

      /**
      * Set all elements of a field to a common scalar: A[i] = s.
      *
      * \param A  field (LHS)
      * \param s  scalar (RHS)
      */
      virtual void eqS(FieldT& A, double s) const = 0;

      /**
      * Compound addition assignment for fields : A[i] += B[i].
      *
      * \param A  field (LHS)
      * \param B  field (RHS)
      */
      virtual void addEqV(FieldT& A, FieldT const & B) const = 0;

      /**
      * Allocate blocks
      */
      virtual void allocateBlocks() = 0; 

   };

   // Public inline member functions

   /*
   * Get derivative of free energy w/ respect to a unit cell parameter.
   */
   template <int D, class PT, class ST>
   inline double MixtureReal<D,PT,ST>::stress(int parameterId) const
   {
      UTIL_CHECK(hasStress_);
      return stress_[parameterId];
   }

   /*
   * Has the stress been computed for the current w fields?
   */
   template <int D, class PT, class ST>
   inline bool MixtureReal<D,PT,ST>::hasStress() const
   {  return hasStress_; }

} // namespace Prdc
} // namespace Pscf
#endif
