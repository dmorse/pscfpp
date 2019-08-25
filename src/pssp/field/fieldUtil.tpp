#ifndef PSSP_FIELD_UTIL_TPP
#define PSSP_FIELD_UTIL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "fieldUtil.h"

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Pssp
{

   template <int D>
   void convertFieldBasisToDft(Basis<D> const & basis,
                               DArray<double> const& components, 
                               RFieldDft<D>& dft)
   {
      // Create Mesh<D> with dimensions of DFT grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star* starPtr; // pointer to current star
      typename Basis<D>::Wave* wavePtr; // pointer to current wave
      std::complex<double> component;   // complex coefficient for star
      std::complex<double> coeff;       // coefficient for wave
      IntVec<D> indices;                // dft grid indices of wave
      int rank;                         // dft grid rank of wave
      int nStar = basis.nStar();        // number of stars
      int is;                           // star index
      int iw;                           // wave index

      is = 0;
      while (is < nStar) {
         starPtr = &(basis.star(is));
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Make real component (coefficient for star basis function)
            component = std::complex<double>(components[is], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis.wave(iw));
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Make complex component for first star
            component = std::complex<double>(components[is], 
                                             -components[is+1]);
            component /= sqrt(2.0);

            // Loop over waves in first star
            starPtr = &(basis.star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis.wave(iw));
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis.star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = conj(component);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis.wave(iw));
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }

            // Increment is by 2 (two stars were processed)
            is += 2;

         } else {
 
            UTIL_THROW("Invalid invertFlag value");
  
         }

      }

   }

   template <int D>
   void convertFieldDfToBasis(Basis<D> const & basis,
                              RFieldDft<D> const & dft, 
                              DArray<double>& components)
   {
      // Create Mesh<D> with dimensions of DFT grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star* starPtr;  // pointer to current star
      typename Basis<D>::Wave* wavePtr;  // pointer to current wave
      std::complex<double> component;    // complex coefficient for star
      IntVec<D> indices;                 // dft grid indices of wave
      int rank;                          // dft grid rank of wave
      int is;                            // star index
      int nStar = basis.nStar();         // number of stars

      // Loop over stars
      is = 0;
      while (is < nStar) {
         starPtr = &(basis.star(is));
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Characteristic wave is first wave of star
            wavePtr = &(basis.wave(starPtr->beginId));
            indices = wavePtr->indicesDft;
            rank = dftMesh.rank(indices);

            // Compute component value
            component = std::complex<double>(dft[rank][0], dft[rank][1]);
            component /= wavePtr->coeff;
            UTIL_CHECK(abs(component.imag()) < 1.0E-8);
            components[is] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either first wave of 1st star or last wave of 2nd star.
            wavePtr = &(basis.wave(starPtr->beginId));
            if (wavePtr->implicit) {
               starPtr = &(basis.star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &(basis.wave(starPtr->endId-1));
               UTIL_CHECK(!(wavePtr->implicit));
            } 
            indices = wavePtr->indicesDft;
            rank = dftMesh.rank(indices);

            // Compute component value
            component = std::complex<double>(dft[rank][0], dft[rank][1]);
            UTIL_CHECK(abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);
            components[is] = component.real();
            components[is+1] = -component.imag();

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }

} // namespace Pscf::Pssp
} // namespace Pscf
#endif
