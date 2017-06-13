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

namespace Pscf{
   template <int D> class Mesh;
   template <int D> class UnitCell;
}

namespace Pscf { 
namespace Pssp
{ 

   using namespace Util;

   /**
   * Basis function for pseudo-spectral scft.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   class Basis {
   
      struct Wave {
         std::complex<double> coeff;
         double sqNorm;
         IntVec<D> vector;
         int starId;
      };
   
      struct Star {
         int size; 
         int beginId; 
         int endId;  
         int invertFlag; 
         int signFlag; 
         IntVec<D> wave;
         bool cancel;
      };

      // Constructor  
      Basis();
 
      // Associate
      void setUnitCell(const UnitCell<D>& unitCell);
      void setMesh(const Mesh<D>& mesh);
   
      // Initialize
      void allocate();
      void deallocate();
   
      // Accessors 
      int nWave() const;
      int nStar() const;
   
      // Get Wave, access by integer index
      Wave& wave(int i);
   
      // Get wave, access wave by IntVec wavevector
      Wave& wave(IntVec<D> vector);
   
      // Get Star, access by integer index
      Star& star(int i);
   
   private:
   
      int nWave_;
      int nStar_;
   
      DArray<Wave> waves_;
      DArray<Star> stars_;
   
      // Indexing that allows identification by IntVec
      DArray<int> waveId_;
      
      UnitCell<D> const * unitCell_;
      Mesh<D> const * mesh_;
   };

}
}
#include "Basis.tpp"
#endif
