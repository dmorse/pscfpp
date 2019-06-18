#ifndef PSSP_BASIS_H
#define PSSP_BASIS_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <pssp/field/RFieldDft.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/containers/DMatrix.h>
#include <set>

namespace Pscf { 

   template <int D> class UnitCell;

namespace Pssp
{ 

   using namespace Util;

   /**
   * Symmetry-adapted basis for pseudo-spectral scft.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   class Basis 
   {

   public:

      /**
      * Wavevector used to construct a basis function.
      */ 
      class Wave 
      {

      public:

         // Coefficient of wave within star basis function
         std::complex<double> coeff;

         // Square magnitude of associated wavevector
         double sqNorm;

         // Integer indices of wave, on discrete Fourier transform mesh.
         IntVec<D> indicesDft;

         // Integer indices of wave, in first Brillouin zone.
         IntVec<D> indicesBz;

         // Index of star containing this wavevector
         int starId;

         // Is this wave represented implicitly in DFT of real field?
         bool implicit;

      };

      /**
      * List of wavevectors that are related by space-group symmetries.
      *
      * The indices of the wavevectors in a star form a continuous block. 
      * Within this block, waves are listed in descending lexigraphical 
      * order of their integer (ijk) indices, with more signficant indices 
      * listed first.
      */ 
      class Star 
      {

      public:

         /**
         * Eigenvalue of negative Laplacian for this star.
         *
         * Equal to square norm of any wavevector in this star.
         */
         double eigen;

         /**
         * Array of derivatives of eigenvalue w/ respect to lattice parameters.
         */
         FArray<double, 6> dEigen;

         /**
         * Number of wavevectors in the star.
         */
         int size; 

         /**
         * Wave index of first wavevector in star.
         */
         int beginId; 

         /**
         * Wave index of last wavevector in star.
         */
         int endId;  

         /**
         * Index for inversion symmetry of star.
         *
         * A star is said to be closed under inversion iff, for each vector
         * G in the star, -G is also in the star. If a star S is not closed 
         * under inversion, then there is another star S' that is related to 
         * S by inversion, i.e., such that for each G in S, -G is in S'. 
         * Stars that are related by inversion are listed consecutively.
         * 
         * If a star is closed under inversion, then invertFlag = 0.
         *
         * If a star is not closed under inversion, then invertFlag = +1 
         * or -1, with inverFlag = +1 for the first star in the pair of 
         * stars related by inversion and invertFlag = -1 for the second.
         *
         * In a centro-symmetric group, all stars are closed under 
         * inversion. In a non-centro-symmetric group, some stars may 
         * still be closed under inversion.
         */
         int invertFlag; 

         /**
         * Integer indices of characteristic wave of this star.
         *
         * Wave given here is in or on boundary of first Brillouin zone.
         * As a result, computing the norm of this wave must yield eigen.
         * For invertFlag = 0 or 1, this is the first wave in the star.
         * For invertFlag = -1, this is the last wave in the star.
         */
         IntVec<D> waveBz;

         /**
         * Is this star cancelled, i.e., associated with a zero function?
         *
         * The cancel flag is true iff there is no nonzero basis function
         * associated with this star.
         */
         bool cancel;

      };

      // Public member functions of Basis

      /**
      * Default constructor.
      */
      Basis();

      /**
      * Construct basis for a specific grid and space group.
      *
      * Proposal: Initially implementation functions correctly only for
      * identity group, withgroupName == 'I'. 
      */
      void makeBasis(const Mesh<D>& mesh, const UnitCell<D>& unitCell, 
                     std::string groupName);

      /**
      * Construct basis for a specific grid and space group.
      */
      void makeBasis(const Mesh<D>& mesh, const UnitCell<D>& unitCell, 
                     const SpaceGroup<D>& group);

      /**
      * Convert field from symmetry-adapted representation to DFT.
      *
      * \param components coefficients of symmetry-adapted basis functions.
      * \param dft complex DFT representation of a field.
      */
      void convertFieldComponentsToDft(DArray<double>& components, 
                                       RFieldDft<D>& dft);   

      /**
      * Convert DFT of real field to symmetry-adapted representation.
      *
      * \param dft complex DFT representation of a field.
      * \param components coefficients of symmetry-adapted basis functions.
      */
      void convertFieldDftToComponents(RFieldDft<D>& dft, 
                                       DArray<double>& components);   

      /**
      * Update values after change in unit cell parameters.
      */
      void update();

      // Accessors

      /**
      * Total number of wavevectors.
      */
      int nWave() const;

      /**
      * Total number of stars.
      */
      int nStar() const;

      /**
      * Total number of nonzero symmetry-adapted basis functions.
      */
      int nBasis() const;

      /** 
      * Get a specific Wave, access by integer array index.
      */
      Wave& wave(int i);

      /** 
      * Get a Wave, access wave by IntVec of indices.
      */
      Wave& wave(IntVec<D> vector);

      /** 
      * Get a Star, access by integer index.
      */
      Star& star(int i);

   private:

      /// Array of all Wave objects (all wavevectors)
      DArray<Wave> waves_;

      /// Array of Star objects (all stars of wavevectors).
      GArray<Star> stars_;

      /// Indexing that allows identification by IntVec
      DArray<int> waveIds_;

      /// Total number of wavevectors
      int nWave_;

      /// Total number of stars.
      int nStar_;

      /// Pointer to associated UnitCell<D>
      const UnitCell<D>* unitCellPtr_;

      /// Pointer to associated Mesh<D>
      const Mesh<D>* meshPtr_;

      /**
      * Wave struct used for internal processing.
      */
      struct NWave {
         double sqNorm;
         IntVec<D> indicesDft;
         IntVec<D> indicesBz;
      };

      struct NWaveComp {
         bool operator() (NWave, NWave) const;
      };

      /**
      * Construct array of ordered waves.
      */
      void makeWaves();

      /**
      * Sort waves of equal magnitude into stars related by symmetry.
      */
      void makeStars(const SpaceGroup<D>& group);

      /**
      * Access associated Mesh<D> as reference.
      */
      const Mesh<D>& mesh() { return *meshPtr_; }

      /**
      * Access associated UnitCell<D> as reference.
      */
      const UnitCell<D>& unitCell() { return *unitCellPtr_; }

   };

   template <int D>
   inline int Basis<D>::nWave() const
   {  return nWave_; }

   template <int D>
   inline int Basis<D>::nStar() const
   {  return nStar_; }

   template <int D>
   inline 
   typename Basis<D>::Wave& Basis<D>::wave(int i)
   {  return waves_[i]; }

   template <int D>
   inline 
   typename Basis<D>::Star& Basis<D>::star(int i)
   {  return stars_[i]; }

} // namespace Pscf:Pssp
} // namespace Pscf
#include "Basis.tpp"
#endif
