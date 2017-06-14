#ifndef PSSP_BASIS_H
#define PSSP_BASIS_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <util/containers/DArray.h>

namespace Pscf { 
namespace Pssp
{ 
   template <int D> class Mesh;
   template <int D> class UnitCell;
   template <int D> class RFieldDft;

   using namespace Util;

   /**
   * Basis function for pseudo-spectral scft.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   class Basis {
  
      /**
      * A wavevector used in the construction of symmetry adapted basis functions.
      */ 
      class Wave {

         // Coefficient of this wave within associated star basis function
         std::complex<double> coeff;

         // Square magnitude of associated wavevector
         double sqNorm;

         // Integer indices of this wavevector
         IntVec<D> indices;

         // Index of star containing this wavevector
         int starId;

         // Is this wave represented implicitly in DFT of real field?
         bool implicit;

      };
  
      /**
      * A list of wavevectors that are related by space-group symmetry operations.
      *
      * The indices of the wavevectors in a star form a continuous block. Within
      * this block, waves are listed in descending lexigraphical order of their 
      * integer (ijk) indices, with more signficant indices listed first.
      */ 
      class Star 
      {

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
         * Index for symmetry of star under inversion.
         *
         * A star is said to be closed under inversion iff, for each vector G
         * in the start, -G is also in the star. If a star S is not closed 
         * under inversion, then there is another star S that is related to 
         * S by inversion, i.e., such that for each G in S, -G is in S'. Stars
         * that are related by inversion are always listed consecutively.
         * 
         * If a star is closed under inversion, invertFlag = 0.
         *
         * If a start is not closed under inversion, then invertFlag = +1 or -1,
         * with inverFlag = +1 for the first star in the pair of stars related 
         * by inversion and invertFlag = -1 for the second.
         *
         * In a centro-symmetric group, all stars are closed under inversion.
         * In a non-centro-symmetric group, some stars may still be closed under
         * inversion.
         */
         int invertFlag; 

         /**
         * Index for symmetry of associated basis function under inversion.
         *
         * If basis function is even under inversion signFlag = 1.
         * If basis function is odd under inversion signFlag = -1.
         */
         int signFlag; 

         /**
         * Integer indices of characteristic wave of this star.
         *
         * For invertFlag = 0 or 1, this is the first wave in the star.
         * For invertFlag = -1, this is the last wave in the star.
         */
         IntVec<D> wave;

         /**
         * Is this star cancelled, i.e., associated with a zero function?
         *
         * The cancel flag is true iff there is not a nonzero basis function
         * associated with this star.
         */
         bool cancel;

      };

      /**
      * Default constructor.
      */
      Basis();
 
      /**
      * Construct basis for a specific grid and space group.
      *
      * Proposal: Initially implementation functions correctly only if groupName == 'I'. 
      */
      void makeBasis(IntVec<D> meshDimensions, const UnitCell<D>& unitCell, std::string groupName);
     
      /**
      * Convert field from symmetry-adapted representation to complex DFT.
      *
      * \param components coefficients of symmetry-adapted basis functions.
      * \param dft complex DFT representation of a field.
      */
      void convertFieldComponentsToDFT(DArray<double> components, RFieldDft<D>& dft);   

      /**
      * Convert DFT of real field to symmetry-adapted representation.
      *
      * \param dft complex DFT representation of a field.
      * \param components coefficients of symmetry-adapted basis functions.
      */
      void convertFieldDftToComponents(RFieldDft<D>& dft, DArray<double> components);   

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
      DArray<Star> stars_;
   
      /// Indexing that allows identification by IntVec
      DArray<int> waveId_;
     
      /// Total number of wavevectors
      int nWave_;

      /// Total number of stars.
      int nStar_;

      /// Pointer to associated UnitCell<D>
      UnitCell<D>* unitCellPtr_;

      /// Dimensions of associated spatial grid.
      IntVec<D> meshDimensions_;
   
   };

}
}
#endif
