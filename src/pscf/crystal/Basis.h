#ifndef PSSP_BASIS_H
#define PSSP_BASIS_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/SpaceGroup.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>

namespace Pscf { 

   using namespace Util;

   /**
   * Symmetry-adapted basis for pseudo-spectral scft.
   *
   * A symmetry-adapted basis for a periodic structure with a specified 
   * space group symmetry may be represented by a set of basis functions
   * that are each associated with a "star" of wavevectors. A star is a
   * set of reciprocal lattice wavevectors that are related by symmetry
   * operations. Wavevectors within a star all have the same magnitude. 
   * For example, in any 3D cubic space group, the wavevectors in each 
   * star have coefficients {ijk} that are related by changes in sign of 
   * individual coefficients and/or permutations of the three indices. 
   * Each basis function in a symmetry-adapted Fourier basis may be 
   * expressed as a linear superposition of complex exponential plane
   * wave phasors associated with different wavevectors in a single star. 
   * The complex coefficients associated with different plane waves that
   * contribute to a basis function are dictated to within an overall
   * phase by requirement that the basis function be invariant under all 
   * symmetry elements of the space group, and by a normalization condition 
   * requiring that the sum of the squares of the coefficients of waves
   * contributing to a basis function or star equal unity. The symmetry
   * condition implies that coefficients of all waves in a star must be
   * complex numbers of equal magnitude, with phase relationships that
   * are dictated by the space group operations. 
   *
   * Individual wavevectors and stars are represented by instances of 
   * the local classes Basis::Wave and Basis::Star, respectively. After 
   * construction of a basis is completed by calling the makeBasis
   * function, a Basis has a private array of waves (i.e., instances 
   * of Basis::Wave) and an array stars (instances of Basis::Star),
   * which we will refer to in this documentation as the wave array
   * and the star array. (The private member variable names are waves_ 
   * and stars_, respectively).  Elements of the wave array may be 
   * accessed by const reference by the function wave(int id) function, 
   * and elements of the star array may be acccess by the function 
   * star(int id).
   *
   * Elements of the wave array are listed in order of non-decreasing
   * wavevector magnitude, with wavevectors in the same star listed as
   * as a consecutive block. Wavevectors within each star are listed 
   * in order of decreasing order as determined by the integer indices,
   * as defined by the member function Wave::indicesBz, with more
   * signficant digits on the left. For example, the waves of the {111} 
   * star of a cubic structure will be listed in the order:
   *
   *     1   1   1  (first)
   *     1   1  -1
   *     1  -1   1
   *     1  -1  -1
   *    -1   1   1  
   *    -1   1  -1
   *    -1  -1   1
   *    -1  -1  -1  (last)
   *
   * Further information about each wavevector or each star is contained
   * public member variables of the Wave::Wave and Wave::Star classes.
   *
   * \ingroup Pscf_Crystal_Module
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

         /**
         * Coefficient of wave within the associated star basis function.
         */
         std::complex<double> coeff;

         /**
         * Square magnitude of associated wavevector
         */
         double sqNorm;

         /**
         * Integer indices of wave, on a discrete Fourier transform mesh.
         *
         * Components of this IntVec<D> are non-negative. Component i lies
         * within the range [0, meshDimension(i)], where meshDimension(i) 
         * is the number of grid points.
         */
         IntVec<D> indicesDft;

         /**
         * Integer indices of wave, in first Brillouin zone.
         *
         * Components of this IntVec<D> may be negative or positive, and
         * differ from corresponding components of indicesDft by integer
         * multiples of a corresponding mesh dimension, so that indicesBz
         * and indicesDft correspond to equivalent aliased wavevectors 
         * for functions evaluated on the chosen discretization mesh. 
         * The shifts relative to indicesDft are chosen so as to minimize 
         * the norm of the Cartesian wavevector constructed by taking a 
         * linear combination of Bravais lattice basis vectors multiplied 
         * by components of indicesBz.
         */
         IntVec<D> indicesBz;

         /**
         * Index of the star that contains this wavevector.
         */
         int starId;

         /**
         * Is this wave represented implicitly in DFT of real field?
         * 
         * In the discrete Fourier transform (DFT) of a real function, 
         * the coefficients of nonzero wavevector G and -G must be 
         * complex conjugates.  As a result, only one of these two 
         * coefficients is stored in the container RFieldDft used to
         * store the DFT of a real function. For each such pair, the
         * member variable Basis::Wave::implicit is set false for the
         * the wavevector whose coefficient is represented explicitly,
         * and is set true for the wavevector whose coefficient is
         * left implicit.
         */
         bool implicit;

      };

      /**
      * A list of wavevectors that are related by space-group symmetries.
      *
      * The wavevectors in a star form a continuous block within the array
      * waves defined by the Basis classes.  Within this block, waves are 
      * listed in descending lexigraphical order of their integer (ijk) 
      * indices as given by Basis::Wave::indexBz.
      */
      class Star 
      {

      public:

         /**
         * Square magnitude of any wavevector in this star.
         *
         * Eigenvalue of negative Laplacian for the basis function 
         * associated with this star.
         */
         double eigen;

         /**
         * Number of wavevectors in this star.
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
         * under inversion, then there is another star S' that is related 
         * to S by inversion, i.e., such that for each G in S, -G is in S'. 
         * Stars that are related by inversion are listed consecutively
         * within the waves array of the Basis class.
         * 
         * If a star is closed under inversion, then invertFlag = 0.
         *
         * If a star is not closed under inversion, then invertFlag = +1 
         * or -1, with inverFlag = +1 for the first star in the pair of 
         * stars related by inversion and invertFlag = -1 for the second.
         *
         * In a centro-symmetric group (i.e., a group that includes the
         * inversion operation, with an inversion center at the origin
         * of the coordinate system), all stars are closed under 
         * inversion. In this case, all stars in the basis must thus
         * have starInvert = 0.
         *
         * In a non-centro-symmetric group, a stars may either be closed
         * closed under inversion (starInvert = 0) or belong to a pair 
         * of consecutive stars that are related by inversion
         * (starInver = +1 or -1).
         */
         int invertFlag; 

         /**
         * Integer indices indexBz of a characteristic wave of this star.
         *
         * Wave given here is the value of indexBz for a characteristic
         * wave for the star.  For invertFlag = 0 or 1, the characteristic
         * wave is the first wave in the star.  For invertFlag = -1, this 
         * is the last wave in the star.
         */
         IntVec<D> waveBz;

         /**
         * Is this star cancelled, i.e., associated with a zero function?
         *
         * The cancel flag is true iff no nonzero basis function can be
         * constructed from the waves of this star that is invariant under
         * all elements of the space group.  Basis functions are thus 
         * associated only with stars for which cancel == false.
         *
         * Each cancelled star contains a set of waves that are related
         * by symmetry for which the X-ray or neutron scattering intensity 
         * would exhibit a systematic cancellation imposed by the space
         * group symmetry. 
         */
         bool cancel;

      };

      // Public member functions

      /**
      * Default constructor.
      */
      Basis();

      /**
      * Destructor.
      */
      ~Basis();

      /**
      * Construct basis for a specific mesh and space group.
      *
      * This function implements the algorithm for constructing a 
      * basis. It is called internally by the overloaded makeBasis 
      * function that takes a file name as an argument rather than a 
      * SpaceGroup<D> object.
      *
      * \param mesh  spatial discretization grid
      * \param unitCell  crystallographic unitCell
      * \param group  crystallographic space group
      */
      void makeBasis(Mesh<D> const & mesh, 
                     UnitCell<D> const & unitCell, 
                     SpaceGroup<D> const & group);

      /**
      * Construct basis for a specific mesh and space group name.
      *
      * This function attempts to identify a file with a name given by
      * groupName that contains a listing of all of the symmetry elements
      * in the space group. It looks first for a file with the specified
      * name defined relative to the current working directory. Failing
      * that, it looks for a file with the specified name in a standard
      * data directory tree that contains space group files for all 
      * possible space groups (i.e,. all 230 3D groups, 17 2D groups and 
      * 2 1d groups) with file names derived from international table 
      * names. Those files are located in the data/groups directory 
      * tree of the package root directory.
      *
      * If a file containing a valid group is found, this function passes 
      * the resulting SpaceGroup<D> object to the overloaded function
      * makeBasis(Mesh<D> const&, UnitCell<D> const&, SpaceGroup<D> const&).
      *
      * This function throws an exception if a file with specified name 
      * is not found or does not contain a valid space group. 
      *
      * \param mesh  spatial discretization grid
      * \param unitCell  crystallographic unitCell
      * \param groupName  string identifier for the space group
      */
      void makeBasis(Mesh<D> const & mesh, 
                     UnitCell<D> const & unitCell, 
                     std::string groupName);

      /**
      * Update values after change in unit cell parameters.
      */
      void update();

      /**
      * Print a list of all waves to an output stream.
      *
      * \param out output stream to which to write
      * \param outputAll output cancelled waves only if true
      */
      void outputWaves(std::ostream& out, bool outputAll = false) const;

      /**
      * Print a list of all stars to an output stream.
      *
      * \param out output stream to which to write
      * \param outputAll output cancelled waves only if true
      */
      void outputStars(std::ostream& out, bool outputAll = false) const;

      /**
      * Returns true if valid, false otherwise.
      */
      bool isValid() const;

      // Accessors

      /**
      * Total number of wavevectors.
      */
      int nWave() const;

      /**
      * Total number of wavevectors in uncancelled stars.
      */
      int nBasisWave() const;

      /**
      * Total number of stars.
      */
      int nStar() const;

      /**
      * Total number of nonzero symmetry-adapted basis functions.
      */
      int nBasis() const;

      /** 
      * Get a Star, access by integer index.
      */
      Star const & star(int i) const;

      /** 
      * Get a specific Wave, access by integer index.
      */
      Wave const & wave(int i) const;

      /** 
      * Get the integer index of a wave, as required by wave(int i).
      * 
      * This function returns the index of a wavevector within the wave  
      * array, as required as in input to function Wave::wave(int). The 
      * components of the input vector are shifted to a value of indexDft 
      * that lies within the DFT mesh before the index is looked up. 
      *
      * \param vector vector of integer indices of a wave vector.
      */
      int waveId(IntVec<D> vector) const;

   private:

      /// Array of all Wave objects (all wavevectors)
      DArray<Wave> waves_;

      /// Array of Star objects (all stars of wavevectors).
      GArray<Star> stars_;

      /// Indexing that allows identification by IntVec
      DArray<int> waveIds_;

      /// Total number of wavevectors
      int nWave_;

      /// Total number of wavevectors in uncancelled stars
      int nBasisWave_;

      /// Total number of stars.
      int nStar_;

      /// Total number of basis functions (or uncancelled stars).
      int nBasis_;

      /// Pointer to associated UnitCell<D>
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D>
      Mesh<D> const * meshPtr_;

      /**
      * Construct array of ordered waves.
      */
      void makeWaves();

      /**
      * Sort waves of equal magnitude into stars related by symmetry.
      */
      void makeStars(const SpaceGroup<D>& group);

      /**
      * Access associated Mesh<D> as const reference.
      */
      Mesh<D> const & mesh() const { return *meshPtr_; }

      /**
      * Access associated UnitCell<D> as const reference.
      */
      UnitCell<D> const & unitCell() const 
      { return *unitCellPtr_; }

   };

   template <int D>
   inline int Basis<D>::nWave() const
   {  return nWave_; }

   template <int D>
   inline int Basis<D>::nBasisWave() const
   {  return nWave_; }

   template <int D>
   inline int Basis<D>::nStar() const
   {  return nStar_; }

   template <int D>
   inline 
   typename Basis<D>::Wave const & Basis<D>::wave(int i) const
   {  return waves_[i]; }

   template <int D>
   inline 
   typename Basis<D>::Star const & Basis<D>::star(int i) const
   {  return stars_[i]; }

   template <int D>
   int Basis<D>::waveId(IntVec<D> vector) const
   {
      meshPtr_->shift(vector);
      int rank = mesh().rank(vector);
      return waveIds_[rank];
   }

   #ifndef PSSP_BASIS_TPP
   extern template class Basis<1>;
   extern template class Basis<2>;
   extern template class Basis<3>;
   #endif

} // namespace Pscf

#endif
