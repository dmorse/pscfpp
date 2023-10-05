#ifndef PSCF_BASIS_H
#define PSCF_BASIS_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>             // inline waveId
#include <pscf/mesh/Mesh.h>               // inline waveId
#include <util/containers/DArray.h>       // member
#include <util/containers/GArray.h>       // member

namespace Pscf { 
namespace Prdc { 

   template <int D> class UnitCell;
   template <int D> class SpaceGroup;

   using namespace Util;

   /**
   * Symmetry-adapted Fourier basis for pseudo-spectral scft.
   *
   * <b> %Basis Functions, Waves, and Stars </b>:
   *
   * This class constructs a symmetry adapted Fourier basis for periodic
   * structures with a specified space group symmetry. Each basis function
   * is such an expansion is a periodic function that is invariant under 
   * translations by any Bravais lattice transtlation and under any 
   * symmetry operations of the specified space group. Each such basis 
   * function is also an eigenfunction of the Laplacian.
   *
   * A Fourier series expansion of a periodic function that is invariant 
   * under under all translations in a specified Bravais lattice is given
   * by expansion of plane waves with wavevectors that all belong to the 
   * corresponding reciprocal lattice.  Components of all wavevectors 
   * are represented within the Basis class using reciprocal lattice basis 
   * vectors as a basis. In this coordinate system, the coordinates of
   * of each reciprocal lattice vector (i.e., any allowed wavevector in 
   * the Fourier expansion of a periodic function) is given by a list of 
   * D integers, where D is the dimension of space. We often refer to 
   * the integer coefficients of a reciprocal lattice vector in this 
   * basis as the indices of the vector.
   *
   * Each basis function in a symmetry-adapated Fourier basis is equal
   * to a linear superposition of complex exponential plan waves 
   * associated with a set of reciprocal lattice vectors that are related 
   * by space group symmetry operations.  We refer to such a set of 
   * symmetry related wavevectors throughout the code and documentation 
   * as a "star".  For example, in any 3D cubic space group, the 
   * wavevectors in each star have coefficients {ijk} that are related 
   * by changes in sign of individual coefficients and/or permutations 
   * of the three indices. Each wave in the Fourier expansion of a 
   * function belongs to one and only one star. The basis function 
   * associated with a star is thus defined by specifying a complex-valued 
   * coefficient for the complex exponential plane wave exp(iG.r) 
   * associated with each wave G in the associated star, where we use 
   * G.r to represent a dot product of wavevector G and position vector r. 
   *
   * <b> %Wave and %Star Arrays: </b>
   *
   * Individual wavevectors and stars are represented by instances of 
   * the local classes Basis::Wave and Basis::Star, respectively. After 
   * construction of a basis is completed by calling the makeBasis 
   * function, a Basis has a private array of waves (i.e., instances of
   * Basis::Wave) and an array stars (instances of Basis::Star), which
   * we will refer to in this documentation as the wave array and the
   * star array. (The actual array containers are private member 
   * variables named waves_ and stars_, respectively). Individual 
   * elements of the wave array may be accessed by const reference by 
   * the accessor function wave(int id) function, and elements of the 
   * star array may be acccess by the function star(int id).
   *
   * Elements of the wave array are listed in order of non-decreasing
   * Cartesian wavevector magnitude, with wavevectors in the same star 
   * listed as as a consecutive block. Wavevectors within each star are 
   * listed in order of decreasing order as determined by the integer 
   * indices, as defined by the member function Wave::indicesBz, with more
   * signficant digits on the left. For example, the waves of the {211} 
   * star of a cubic structure with full cubic symmetry will be listed 
   * in the order:
   * \code
   *     1   1   1  
   *     1   1  -1
   *     1  -1   1
   *     1  -1  -1
   *    -1   1   1  
   *    -1   1  -1
   *    -1  -1   1
   *    -1  -1  -1  
   * \endcode
   * Further information about each wavevector and each star is 
   * contained in public member variables of the Wave::Wave and 
   * Wave::Star classes. See the documentation of these local classes
   * for some further details. 
   *
   * <b> Discrete Fourier Transforms and Wavevector Aliasing </b>:
   *
   * This class is designed for contexts in which the values of real
   * periodic functions are represented numerically by values evaluated 
   * on a regular grid of points within each primitive unit cell of 
   * the crystal, and in which a discrete Fourier transform (DFT) is 
   * used to transform between Fourier and real-space representations 
   * of such function. Construction of a symmetry adapted Fourier 
   * basis requires a description of this spatial discretization grid,
   * given by an instance of class Mesh<D>, in addition UnitCell<D>
   * and SpaceGroup<D> objects. 
   * 
   * Interpretation of indices for wavevectors in this context is
   * complicated by the phenomena of ``aliasing" of wavevectors used in
   * the Fourier expansion of periodic functions that are defined only
   * at the nodes of a regular grid.  In what follows, let N_i denote 
   * the the number of lattice grid points along direction i in a 
   * coordinate system in which Bravais lattice vectors are used as
   * basis vectors to construct corresponding Cartesian vectors. 
   * Any two reciprocal lattice vectors for which the integer indices
   * associate with any direction i differs by an integer multiple of 
   * N_i are equivalent in the sense that the equivalent that they 
   * yield equivalent values for the complex exponential exp(iG.r) 
   * for all values of r that lie on lattice nodes. Such vectors 
   * are referred to as aliases of one another. Two different schemes
   * are used in here to assign a choose a unique choice for a list 
   * of indices for each distinct wavevector. 
   *
   * (1) The "DFT" (discrete Fourier transform) indices of a wavevector 
   * is the choice of a list of indices such that the index associated 
   * with direction i, denoted by m_i, is in the range 0 <= m_i < N_i.
   *
   * (2) The "BZ" (Brillouin zone) indices of a wavevector is a list of
   * indices of the wavevector or one of its aliases chosen such that
   * the the Cartesian norm of the wavevector constructed as a linear
   * superposition of reciprocal lattice vectors multiplied these
   * indices is less than or equal to the norm of any other alias of
   * the wavevector. In cases in which a set of two or more aliases 
   * of a wavevector have the same Cartesian norm, the BZ indices are 
   * chosen to be those of the member of the set for which the indices
   * are "largest" when lists of indices are compared by treating
   * earlier indices are more signficant, as done in the ordering of 
   * waves within a star in the waves array. 
   *
   * The number of distinct waves in the waves array is equal to the 
   * number of grid points in the corresponding spatial mesh. Each wave
   * has a unique list of DFT indices and also a unique list of BZ 
   * indices, which are listed in the "indicesDft" and "indicesBz"
   * IntVec<D> members of the Basis<D>::Wave class. Values of the 
   * square magnitude of wavevectors are always computed using the
   * indices of each wavevector.
   * 
   * <b> Space Group Symmetry </b>:
   *
   * A space group symmetry is defined by a pair (R,t) in which R is a
   * linear transformation that may represented as a matrix and t is a
   * translation vector. The translation vector may be either zero or a
   * vector whose components in a basis of Bravais lattice vectors are
   * all fractions (i.e., translations by fractions of a unit cell).
   * The effect of such a symmetry operation on a position vector r is 
   * to transform r to a tranformed position r' given by r' = R.r + t, 
   * where R.r represents the result of applying linear transformation
   * (or matrix) R to vector r.
   *
   * Applying such a symmetry operation to a plane wave f(r) = exp(iG.r) 
   * with a wavevector G transforms f(r) into a different plane wave 
   * f'(r) given by f'(r) = exp[i G.(Rr + t)] =  exp[i(G'.r + theta)] 
   * with a transformed wavevector G' given by G' = G.R, with a phase
   * shift theta given by the dot product theta = G.t.  This operation
   * can be expressed using matrix notation by expressing all position
   * vectors such as r, r' and t as column vectors, all wavevectors 
   * such as G and G' as row vectors, and a linear transformation 
   * denoted by R by a matrix. When applying symmetry operations to 
   * wavevectors defined using a discrete Fourier transform, the equation
   * G' = G.R is interpreted using the representation of wavevectors 
   * G and G' using their BZ (Brillouin zone) indices.
   *
   * The plane waves in a star all have the same Cartesian magnitude, 
   * as a result of the fact that the linear transformation represented 
   * by R for a symmetry operation (R,t) must be a norm preserving
   * (i.e., unitary) transformation such as a rotation, reflection or
   * inversion. For any such operation, |G| = |G.R|, where |...| 
   * represents a Cartesian norm (i.e., a Euclidean norm defined using 
   * Cartesian representations of G and R).
   *
   * The requirement that a basis function be invariant under all of the 
   * symmetry operations in a space group imposes a set of relationships 
   * among the coefficients of different waves in the expansion of the
   * basis function associated with a star. Suppose G and G' are two
   * wavevectors in the same star that are related by a transformation 
   * or matrix R associated with a symmetry operation S=(R,t), such that 
   * G'=G.R. Let c and c' be the complex coefficients of planes waves 
   * associated with G and G' within the sum that defines a basis 
   * function. The requirement that the basis function be invariant 
   * under the action of symmetry S requires that c' = c exp(iG.t). The 
   * resulting set of requirements among all the waves in a star implies 
   * that the coefficients associated with different waves in a star must 
   * all have the same absolute magnitude, because |exp(iG.t)|=1, and
   * imposes a set of relationships among the phases (complex arguments)
   * of these coefficients. The magnitude of the coefficients of all waves 
   * in a star is uniquely determined in this class by imposing a 
   * normalization condition requiring that the sum of squares of absolute 
   * magnitudes of coefficients of waves in a star be equal to unity, or 
   * that the square magnitude of each coefficient be the inverse of the 
   * number of distinct distinct waves in the star.  The combination of 
   * this normalization condition and the phase conditions imposed by the 
   * requirement of invariance under space group symmetries defines the 
   * basis function associated with each star to within an overall 
   * complex prefactor of unit absolute magnitude (i.e., a phasor).
   *
   * Space group symmetries are represented internally in the class
   * SpaceSymmetry<D> class using a coordinate system in which all position
   * and translation vectors are expanded using the basis of Bravais basis 
   * vectors defined in UnitCell<D>, while wavevectors are expanded using
   * the corresponding reciprocal lattice basis vectors. In this coordinate 
   * system, for any symmetry operation S = (R,t), elements of the matrix
   * representation of the linear transformation R are all integers 
   * (usually 0, +1, or -1) and elements of the translation vector t are 
   * all rational numbers with small denominators (e.g., 1/4, 1/3, 1/2, 
   * etc.) in the range [0,1].
   *
   * <b> Cancelled Stars </b>
   *   
   * A star is said to be "cancelled" if there is no way to construct
   * nonzero basis function from the waves in that star that is invariant
   * under all of the elements of the space group.  The notion of 
   * cancellation is equivalent to the notion of systematic cancellation 
   * used in analysis of Bragg scattering from crystals - waves belonging 
   * to cancelled stars exhibit systematic cancellation in scattering, 
   * and do not yield Bragg peaks. 
   *
   * A star is cancelled iff, for any wavevector G in the star, there 
   * exist two symmetry operations S = (R,t) and S' = (R',t') in the 
   * space group for which G.R = G.R' but G.t != G.t'. It can shown 
   * that if this is satisfied for any wavevector in the star, it must
   * be satisfied for all wavevector in the star. 
   *
   * The Wave class has a bool member variable named cancel that is 
   * set true for cancelled stars and false otherwise. Each basis
   * function in the symmetry-adapted Fourier basis is associated
   * with a non-cancelled star.
   * 
   * The class Star has two integer members named starId and basisId
   * that correspond to different ways of indexing stars. The index
   * starId is the index of the a star within an ordered list of all 
   * stars, including cancelled stars. The starId is also the array 
   * element index of the star with the stars_ array, which includes 
   * cancelled stars. The basisId of an uncancelled star is its index 
   * within a contracted list that contains only uncancelled stars, 
   * skipping over cancelled stars. The basisId is not defined for
   * a cancelled star, and the member variable basisId is set to -1
   * all cancelled stars by convention. A single Star object may be
   * accessed by its starId using the function Basis::star(int starId), 
   * and an uncancelled star may be accessed by its basisId using the 
   * function Basis::basisFunction(int basisId).
   * 
   * <b> Open and Closed Stars </b>
   *
   * Every star is either "open" or "closed" under the action of the
   * inversion operation. A star is closed under inversion if, for 
   * every wavevector G in the star, the negation -G is in the same
   * star.  The negation of a wavevector is defined by inverting all
   * of its BZ indices. A star that is not closed under inversion,
   * or "closed", is "open". 
   *
   * In a basis for a crystal contains an inversion center, for which 
   * which the space group includes a symmetry operation r -> -r + t,
   * all stars are closed. Open stars only apear in space groups that 
   * do not contain an inversion center. A crystal with an inversion
   * center is said to be "centrosymmetric".  Open stars thus only
   * exist in space groups for non-centrosymmetric crystals.
   *
   * Open stars come in pairs that are related by inversion, such that
   * for every wavevector G in one star in the pair, the negation -G is
   * in the other star. Pairs of open stars that are related by inversion
   * are listed consecutively in the Stars array, and correspond to 
   * neighboring blocks of waves in the waves array.
   *
   * The local Star class has an integer member "starInvert" that is
   * set to if a star is closed, 1 if it is the first member of a pair
   * of open stars that are related by inverse, and -1 if the is the 
   * second member of such a pair of stars.
   *
   * Each star is assigned a unique characteristic wave that can be
   * used as a unique identifier for the star.  Local class Star has an
   * IntVec<D> member named waveBz that lists the BZ indices of the
   * characteristic wave of the star. The characteristic wave is taken 
   * to be the first wave in the star for stars with starInvert = 0 or 
   * starInvert = 1. In stars with startInvert = -1, the characteristic
   * wave is taken to be the negation of the characteristic wave of 
   * its partner, which is the previous star.  
   *
   * <b> Phase Conventions for %Basis Functions </b>
   *
   * Closed stars:
   * Phases of the coefficients of waves in each closed star are
   * chosen so that the associated basis function is a real function
   * of position. This is done by choosing coefficients such that 
   * the coefficient associated with each wavevector and its
   * negation are complex conjugates. This requirement, together 
   * with the normalization condition, determines the basis function 
   * to within an overall sign. The sign of the basis function is
   * fixed by requiring that the coefficient of the characteristic
   * wave of the star (the first wave listed) has a non-negative
   * real part, and a negative imaginary part in the special case
   * in which the coefficient must be chosen to be pure imaginary.
   * 
   * Basis functions associated with pairs of stars that are related 
   * by inversion are defined so as to be complex conjugates of one
   * another. This is done by requiring that the coefficient of the
   * characteristic wave of each such star is real and positive.
   *
   * <b> Fourier expansion of a real Function </b>
   *
   * The expansion of a real function with a specified space group
   * symmetry is specified by specifying a real coefficient for each
   * closed star and a pair of real coefficients for each pair of closed
   * stars. The number of required real coefficients is thus the same as
   * the same number as the number of uncancelled stars or basis 
   * functions. By convention, these coefficients should be stored in 
   * an array in which the index of the coefficient of a the basis 
   * function associated with a closed star is the same as the index of
   * the closed star, and the indices of the two real coefficients of 
   * each pair of basis functions arising from open stars that are
   * related by inversion correspond to the indices of these two stars.
   *
   * Suppose two open stars stars that are related by inversion appear
   * in the stars array with indices j and j+1. Let f(r) denote the 
   * basis function associated with star j and let its complex conjugate 
   * f^{*}(r) denote the basis function associated with star j+1. Let
   * A denote an array of real coefficients used to expand a real 
   * periodic basis function. Let a = A[j] and b = A[j+1] denote the 
   * coefficients that appear in corresponding elements j and j + 1
   * of array A. The contribution of stars j and j+1 to the desired
   * real periodic function is given by an expression
   *
   *        (a - ib) f(r) + (a + ib) f^{*}(r) 
   *
   * in which i denotes the square root of -1.
   *
   * The contribution of a basis function f(r) associated with a
   * closed star is simply given by the product of the associated
   * real coefficient and this real basis function. 
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

         /*
         * Default constructor.   
         */
         Wave()
          : coeff(0.0),
            indicesDft(0),
            indicesBz(0),
            starId(0),
            implicit(false),
            sqNorm(0.0)
         {}

      private:

         /**
         * Square magnitude of associated wavevector
         */
         double sqNorm;

         friend class Basis<D>;

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
         * is the negation of the waveBz of the partner star for which
         * invertFlag = 1.
         *
         * The coefficient of wave waveBz is always real and positive.
         */
         IntVec<D> waveBz;

         /**
         * Index of this star in ordered array of all stars.
         *
         * This is the index of this star within array stars_.
         */
         int starId;

         /**
         * Index of basis function associated with this star.
         *
         * If this star is not cancelled (cancel == false), then basisId
         * is set equal to the index of the associated basis function in
         * in a list of nonzero basis functions or (equivalently) in a
         * list of non-cancelled stars.
         * 
         * If this star is cancelled (cancel == true), then basisId = -1.
         */
         int basisId;

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
         * would exhibit a systematic cancellation required by the space
         * group symmetry. 
         */
         bool cancel;

         /*
         * Default constructor.   
         */
         Star()
          : size(0),
            beginId(0),
            endId(0),
            invertFlag(0),
            waveBz(0),
            starId(0),
            basisId(0),
            cancel(false)
         {}

      };

      // Public member functions of Basis<D>

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

      /**
      * Returns true iff this basis is fully initialized.
      */
      bool isInitialized() const;

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
      * Get a specific Wave, access by integer index.
      *
      * \param id  index for a wave vector.
      */
      Wave const & wave(int id) const;

      /** 
      * Get a Star, accessed by integer star index.
      *
      * This function return both cancelled and un-cancelled stars.
      *
      * \param id  index for a star
      */
      Star const & star(int id) const;

      /** 
      * Get an uncancelled Star, accessed by basis function index.
      *
      * \param id  index for a basis function, or an uncancelled star.
      */
      Star const & basisFunction(int id) const;

      /** 
      * Get the integer index of a wave, as required by wave(int id).
      * 
      * This function returns the index of a wavevector within the wave  
      * array, as required as in input to function Wave::wave(int id). 
      * The vector of input components is shifted internally to an 
      * equivalent value that lies within the conventional DFT mesh,
      * in which all components are non-negative integers.
      *
      * \param vector vector of integer indices of a wave vector.
      */
      int waveId(IntVec<D> vector) const;

   private:

      /**
      * Array of all Wave objects (all wavevectors). 
      *
      * Waves are ordered in non-decreasing order by wavevector norm, with
      * waves in the same star order in a consecutive block, with waves in
      * each star ordered by integer indices.
      *
      * The capacity of array stars_ is equal to nWave_, which is also equal
      * to the number of points in the associated spatial mesh.
      */
      DArray<Wave> waves_;

      /**
      * Array of Star objects (all stars of wavevectors).
      * 
      * The final size of array stars_ is equal to nStar_.
      */
      GArray<Star> stars_;

      /**
      * Look-up table for identification of waves by IntVec
      *
      * After allocation, the capacity_ of waveIds_ is equal to nWave_.
      * The array index of each element of waveIds_ corresponds to the
      * rank of the wavevector within a k-space DFT mesh, and the value
      * is the index of the corresponding Wave within the waves_ array.
      */
      DArray<int> waveIds_;

      /**
      * Look-up table for uncancelled stars index by basis function id.
      *
      * After allocation, the capacity_ of starIds_ is equal to nBasis_.
      * The array index of each element of starIds_ corresponds to the
      * index of a basis function, while the element value is the index
      * of the corresponding un-cancelled Star in the stars_ array.
      */
      DArray<int> starIds_;

      /**
      * Total number of wavevectors, including those in cancelled stars.
      *
      * This must equal the number of grid points in the mesh.
      */
      int nWave_;

      /**
      * Total number of wavevectors in uncancelled stars
      */ 
      int nBasisWave_;

      /**
      * Total number of stars, including cancelled stars.
      */
      int nStar_;

      /**
      * Total number of basis functions, or uncancelled stars.
      */
      int nBasis_;

      /**
      * Pointer to an associated UnitCell<D>
      */
      UnitCell<D> const * unitCellPtr_;

      /**
      * Pointer to an associated Mesh<D>
      */
      Mesh<D> const * meshPtr_;

      /**
      * Has this basis been fully initialized?
      */
      bool isInitialized_;

      /**
      * Construct an array of ordered waves.
      */
      void makeWaves();

      /**
      * Identify stars of wavevectors related by symmetry.
      */
      void makeStars(const SpaceGroup<D>& group);

      /**
      * Access associated Mesh<D> as const reference.
      */
      Mesh<D> const & mesh() const 
      { return *meshPtr_; }

      /**
      * Access associated UnitCell<D> as const reference.
      */
      UnitCell<D> const & unitCell() const 
      { return *unitCellPtr_; }

   };

   // Inline functions

   template <int D>
   inline int Basis<D>::nWave() const
   {  return nWave_; }

   template <int D>
   inline int Basis<D>::nBasisWave() const
   {  return nBasisWave_; }

   template <int D>
   inline int Basis<D>::nStar() const
   {  return nStar_; }

   template <int D>
   inline 
   typename Basis<D>::Wave const & Basis<D>::wave(int id) const
   {  return waves_[id]; }

   template <int D>
   inline 
   typename Basis<D>::Star const & Basis<D>::star(int id) const
   {  return stars_[id]; }

   template <int D>
   inline 
   typename Basis<D>::Star const & Basis<D>::basisFunction(int id) const
   {  return stars_[starIds_[id]]; }

   template <int D>
   int Basis<D>::waveId(IntVec<D> vector) const
   {
      meshPtr_->shift(vector);
      int rank = mesh().rank(vector);
      return waveIds_[rank];
   }

   template <int D>
   inline bool Basis<D>::isInitialized() const
   {  return isInitialized_; }

   #ifndef PSCF_BASIS_TPP
   extern template class Basis<1>;
   extern template class Basis<2>;
   extern template class Basis<3>;
   #endif

} // namespace Prdc
} // namespace Pscf
#endif
