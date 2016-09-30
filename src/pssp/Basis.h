#ifndef PSSP_BASIS_H
#define PSSP_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf { 
namespace Pssp
{ 

   using namespace Util;

   /**
   * Basis function for pseudo-spectral scft.
   *
   * \ingroup Pscf_Pssp_Module
   */
   template <int D>
   class Basis {
   
      struc Wave {
         std::complex<double> coeff;
         double sqNorm;
         IntVec<D> vector;
         int starId;
      }
   
      struc Star {
         int size; 
         int beginId; 
         int endId;  
         int invertFlag; 
         int signFlag; 
         IntVec<D> wave;
         bool cancel;
      };

      // Constructor  
      Grid();
 
      // Associate
      void setUnitCell(const UnitCell& unitCell)
      void setGrid(const Grid& grid)
   
      // Initialize
      void allocate();
      void deallocate();
   
      // Accessors 
      int nWave() const;
      int nStar() const;
   
      // Get Wave, access by integer index
      Wave& wave(int i);
   
      // Get wave, access wave by integer vector
      Wave& wave(IntVec<D> vector);
   
      // Get Star, access by integer index
      Star& star(int i);
   
   private:
   
      int nWave_;
      int nStar_;
   
      DArray<Wave> waves_;
      DArray<Wave> stars_;
   
      // Indexing that allows identification by IntVec
      DArray<int> waveId_;   
      
   };

}
}
#endif
