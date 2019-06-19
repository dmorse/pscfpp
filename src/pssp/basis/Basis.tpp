#ifndef PSSP_BASIS_TPP
#define PSSP_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"
#include "TWave.h"
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <vector>
#include <set>

namespace Pscf {
namespace Pssp
{

   /*
   * Constructor.
   */
   template <int D>
   Basis<D>::Basis()
    : nWave_(0), 
      nStar_(0), 
      unitCellPtr_(0), 
      meshPtr_(0)
   {}

   /*
   * Construct basis for pseudo-spectral scft.
   */
   template <int D>
   void Basis<D>::makeBasis(const Mesh<D>& mesh, 
                            const UnitCell<D>& unitCell,
                            std::string groupName)
   {
      SpaceGroup<D> group;
      if (groupName == "I") {
         // Create identity group
         group.makeCompleteGroup();
       } else {
         UTIL_THROW("Unimplemented space group");
       }
       makeBasis(mesh, unitCell, group);
   }

   /*
   * Construct basis for pseudo-spectral scft.
   */
   template <int D>
   void Basis<D>::makeBasis(const Mesh<D>& mesh, 
                            const UnitCell<D>& unitCell,
                            const SpaceGroup<D>& group)
   {
      // Save pointers to mesh and unit cell
      meshPtr_ = &mesh;
      unitCellPtr_ = &unitCell;

      // Allocate arrays
      nWave_ = mesh.size();
      waves_.allocate(nWave_);
      waveIds_.allocate(nWave_); 

      // Make sorted array of waves
      makeWaves();

      // Identify stars of waves that are related by symmetry
      makeStars(group);

   }

   /*
   * Construct ordered list of waves.
   * 
   * On exit:
   *  - Array waves_ contains list of waves ordered by sqNorm.
   *  - Each wave has indicesDft, indicesBz and sqNorm set.
   *  - stars_ array is still empty.
   */
   template <int D>
   void Basis<D>::makeWaves()
   {
      IntVec<D> meshDimensions = mesh().dimensions();
      std::vector< TWave<D> > twaves;
      twaves.reserve(nWave_);

      // Loop over dft mesh to generate all waves, add to twaves
      TWave<D> w;
      IntVec<D> v;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         w.indicesDft = itr.position();
         v = shiftToMinimum(w.indicesDft, meshDimensions, *unitCellPtr_);
         w.indicesBz = v;
         w.sqNorm = unitCell().ksq(v);
         twaves.push_back(w);
      }

      // Sort temporary array twaves
      TWaveNormComp<D> comp;
      std::sort(twaves.begin(), twaves.end(), comp);

      // Copy temporary array twaves into member variable waves_
      for (int i = 0; i < nWave_; ++i) {
         waves_[i].sqNorm = twaves[i].sqNorm;
         waves_[i].indicesDft = twaves[i].indicesDft;
         waves_[i].indicesBz = twaves[i].indicesBz;
      }

   }

   template <int D>
   void Basis<D>::makeStars(const SpaceGroup<D>& group)
   {
      /* 
      * Local containers that hold TWave<D> objects:
      * list : set of waves of equal norm, compared using indicesDft
      * star : set of symmetry-related waves, compared by indicesDft
      * tempStar : temporary star, sorted by descending indicesBz
      * tempList:= temporary list, ordered by star
      */
      std::set< TWave<D>, TWaveDftComp<D> > list;  
      std::set< TWave<D>, TWaveDftComp<D> > star;  
      std::vector< TWave<D> > tempStar;  
      GArray< TWave<D> > tempList;            

      typename std::set< TWave<D>, TWaveDftComp<D> >::iterator rootItr;
      typename std::set< TWave<D>, TWaveDftComp<D> >::iterator setItr;
      TWave<D> wave;
      Basis<D>::Star newStar;

      std::complex<double> coeff;
      double Gsq;
      double Gsq_max = 1.0;
      double epsilon = 1.0E-8;
      double twoPi = 2.0*Constants::Pi;
      IntVec<D> rootVec;
      IntVec<D> vec;
      int listId = 0;      // id for this list
      int listBegin = 0;   // id of first wave in this list
      int listEnd = 0;     // (id of last wave in this list) + 1
      int starId = 0;      // id for this star
      int starBegin = 0;   // id of first wave in this star
      int starSize = 0;    // size (number of waves) in this star
      int i, j, k;
      bool cancel;

      for (i = 0; i <= nWave_; ++i) {

         // Determine if this wave begins a new list
         bool newList = false;
         if (i == nWave_) {
            listEnd = i;
            newList = true;
         } else {
            Gsq = waves_[i].sqNorm;
            if (Gsq > Gsq_max + epsilon) {
               Gsq_max = Gsq;
               listEnd = i;
               newList = true;
            }
         }

         // Process new list
         if (newList && listEnd > 0) {

            // Copy waves of equal norm into container "list"
            list.clear();
            tempList.clear();
            for (j = listBegin; j < listEnd; ++j) {
               wave.indicesDft = waves_[j].indicesDft;
               wave.indicesBz = waves_[j].indicesBz;
               wave.sqNorm = waves_[j].sqNorm;
               list.insert(wave);
            }

            // Loop over stars within this list
            IntVec<D> nVec;
            rootItr = list.begin();
            int nextInvert = 1;
            while (list.size() > 0) {

               // Construct a star from root vector.
               rootVec = rootItr->indicesBz;
               Gsq = unitCell().ksq(rootVec);
               cancel = false;
               star.clear();
               for (j = 0; j < group.size(); ++j) {

                  // Apply symmetry (i.e., multiply by rotation matrix)
                  vec = rootVec*group[j];

                  UTIL_CHECK(abs(Gsq - unitCell().ksq(vec)) < epsilon);

                  // Initialize TWave object
                  wave.sqNorm = Gsq;
                  wave.indicesBz = vec;
                  wave.indicesDft = vec;
                  mesh().shift(wave.indicesDft);
                  wave.phase = 0.0;

                  // Compute phase for basis function coefficient
                  wave.phase = 0.0;
                  for (k = 0; k < D; ++k) {
                     wave.phase += wave.indicesBz[k]*(group[j].t(k));
                  }
                  while (wave.phase > 0.5) {
                     wave.phase -= 1.0;
                  }
                  while (wave.phase <= -0.5) {
                     wave.phase += 1.0;
                  }
                  wave.phase *= twoPi;

                  // Check for cancellation of star
                  if (vec == rootVec) {
                     if (abs(wave.phase) > 1.0E-6) {
                        cancel = true;
                        wave.phase = 0.0;
                     }
                  }

                  // Add wave to star (insertion fails if not unique)
                  star.insert(wave);
               }
               starSize = star.size();

               // Append waves in star to tempStar, erase from list
               tempStar.clear();
               setItr = star.begin(); 
               for ( ; setItr != star.end(); ++setItr) {
                  list.erase(*setItr);
                  tempStar.push_back(*setItr);
               }

               // Sort tempStar vector, ordered by indicesBz
               TWaveBzComp<D> waveBzComp;
               std::sort(tempStar.begin(), tempStar.end(), waveBzComp);
         
               // Append contents of ordered star tempStar to tempList 
               for (j = 0; j < tempStar.size(); ++j) {
                  tempList.append(tempStar[j]);
               }
               UTIL_CHECK(tempList.size()+list.size() == listEnd-listBegin);

               // Initialize new Star object
               newStar.eigen = Gsq;
               newStar.beginId = starBegin;
               newStar.endId = newStar.beginId + star.size();
               newStar.size = star.size();
               newStar.cancel = cancel;

               // Determine newStar.invertFlag, find root of next star
               if (nextInvert == -1) {

                  // If this star is the 2nd of pair related by symmetry,
                  // set root of next star to first wave of remaining list.

                  newStar.invertFlag = -1;
                  rootItr = list.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the second of a pair of partners,
                  // then determine if it is closed under inversion.

                  // Compute negation of root vector, shift to DFT mesh
                  nVec.negate(rootItr->indicesBz);
                  (*meshPtr_).shift(nVec);
       
                  // Search for negation of root vector within star.
                  bool negationFound = false;
                  setItr = star.begin(); 
                  for ( ; setItr != star.end(); ++setItr) {
                     if (nVec == setItr->indicesDft) {
                        negationFound = true;
                        break;
                     }
                  }

                  if (negationFound) {

                     // If star is closed under inversion, then root
                     // of next star is first vector of remaining list.

                     newStar.invertFlag = 0;
                     rootItr = list.begin();
                     nextInvert = 1;

                  } else {

                     newStar.invertFlag = 1;
                     nextInvert = -1;

                     // If star is not closed, find negation of root in
                     // the remaining list, and use this negation as root 
                     // of the next star.

                     setItr = list.begin(); 
                     for ( ; setItr != list.end(); ++setItr) {
                        if (nVec == setItr->indicesDft) {
                           negationFound = true;
                           rootItr = setItr;
                           break;
                        }
                     }
                     // On exit after break, rootItr = &nVec

                     if (!negationFound) {
                        UTIL_THROW("Negative of root vector not found");
                     }
                  }

               }

               stars_.append(newStar);
               ++starId;
               starBegin = newStar.endId;
            }
            UTIL_CHECK(list.size() == 0);
            UTIL_CHECK(tempList.size() == listEnd - listBegin);

            // Copy tempList array into corresponding section of waves_
            for (j = 0; j < tempList.size(); ++j) {
               k = j + listBegin;
               waves_[k].indicesDft = tempList[j].indicesDft;
               waves_[k].indicesBz = tempList[j].indicesBz;
               waves_[k].sqNorm = tempList[j].sqNorm;
               coeff = std::complex<double>(0.0, tempList[j].phase);
               waves_[k].coeff = exp(coeff);
            }
            // At this point, waves_[k].coeff has unit absolute magnitude.
            // Coefficients are normalized in final processing of stars.

            ++listId;
            listBegin = listEnd;
         }
      }
      nStar_ = stars_.size();

      // Final processing of stars
      double snorm;
      int waveId;
      for (i = 0; i < nStar_; ++i) {

         // Set characteristic wavevector waveBz for each star
         if (stars_[i].invertFlag == -1) {
           waveId = stars_[i].endId - 1;
         } else {
           waveId = stars_[i].beginId;
         }
         stars_[i].waveBz = waves_[waveId].indicesBz;

         // Set starId for all associated waves
         for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
            waves_[j].starId = i;
         }

         // Normalize coeff for all associated waves
         if (!stars_[i].cancel) {
            snorm = 1.0/sqrt(double(stars_[i].size));
            for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
               waves_[j].coeff *= snorm;
            }
         }

         // Set elements of dEigen
         {
            vec = stars_[i].waveBz;
            double element, dEigen;
            int p, q;
            for (j = 0; j < unitCell().nParams(); ++j) {
               dEigen = 0;
               for (p = 0; p < D; ++p){
                  for (q = 0; q < D; ++q){
                     element = unitCellPtr_->dkkBasis(j, p, q);
                     dEigen += vec[p]*vec[q]*element;
                  }
               }
               stars_[i].dEigen[j] = dEigen;
            }
         }

      }

      // Final processing of waves
      IntVec<D> meshDimensions = mesh().dimensions();
      for (i = 0; i < nWave_; ++i) {
         vec = waves_[i].indicesDft;

         // Validity check 
         for (j = 0; j < D; ++j) {
            UTIL_CHECK(vec[j] >= 0);
            UTIL_CHECK(vec[j] < meshDimensions[j]);
         }

         // Set implicit attribute
         if ((vec[D-1] + 1) > (meshDimensions[D-1]/2 + 1)) {
            waves_[i].implicit = true;
         } else {
            waves_[i].implicit = false;
         }

         // Look up table for waves
         waveIds_[mesh().rank(vec)] = i;
      }
  
      #if 0
      // Output all waves
      std::cout << std::endl;
      std::cout << "Waves:" << std::endl;
      for (i = 0; i < nWave_; ++i) {
         std::cout << Int(i,4);
         std::cout << Int(waves_[i].starId, 4);
         std::cout << " |";
         for (j = 0; j < D; ++j) {
            std::cout << Int(waves_[i].indicesDft[j], 4);
         }
         std::cout << " |";
         for (j = 0; j < D; ++j) {
            std::cout << Int(waves_[i].indicesBz[j], 4);
         }
         std::cout << " | ";
         std::cout << Dbl(waves_[i].sqNorm, 12);
         std::cout << "  " << Dbl(waves_[i].coeff.real(), 10);
         std::cout << "  " << Dbl(waves_[i].coeff.imag(), 10);
         std::cout << std::endl;
      }
      #endif
 
      #if 0 
      // Output all stars
      std::cout << std::endl;
      std::cout << "Stars:" << std::endl;
      for (i = 0; i < nStar_; ++i) {
         std::cout << Int(i, 4)
                   << Int(stars_[i].size, 3)
                   << Int(stars_[i].beginId, 5)
                   << Int(stars_[i].endId, 5)
                   << Int(stars_[i].invertFlag, 3)
                   << Int(stars_[i].cancel, 3);
         std::cout << " |";
         for (j = 0; j < D; ++j) {
            std::cout << Int(stars_[i].waveBz[j], 4);
         }
         std::cout << " | " << Dbl(stars_[i].eigen, 12);
         std::cout << std::endl;
      }
      #endif
  
   }

   template <int D>
   void Basis<D>::update()
   {
      IntVec<D> vec;
      int i, j, p, q;
      double element, dEigen;

      // Process stars
      for (i = 0; i < nStar_; ++i) {
         vec = stars_[i].waveBz;
         stars_[i].eigen = unitCell().ksq(vec);

         for (j = 0; j < unitCell().nParams(); ++j) {
            dEigen = 0.0;
            for (p = 0; p < D; ++p){
               for (q = 0; q < D; ++q){
                  element = unitCell().dkkBasis(j, p, q);
                  dEigen += vec[p]*vec[q]*element;
               }
            }
            stars_[i].dEigen[j] = dEigen;
         }

      }

      // Process waves
      for (i = 0; i < nWave_; ++i) {
         vec = waves_[i].indicesBz;
         waves_[i].sqNorm = unitCell().ksq(vec);
      }

   }

   template <int D>
   void Basis<D>::convertFieldComponentsToDft(DArray<double>& components, 
                                              RFieldDft<D>& dft)
   {
      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      Star* starPtr;                  // pointer to current star
      Wave* wavePtr;                  // pointer to current wave
      std::complex<double> component; // coefficient for star
      std::complex<double> coeff;     // coefficient for wave
      IntVec<D> indices;              // dft grid indices of wave
      int rank;                       // dft grid rank of wave
      int is;                         // star index
      int iw;                         // wave index

      is = 0;
      while (is < nStar_) {
         starPtr = &stars_[is];
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Make real component (coefficient for star basis function)
            component = std::complex<double>(components[is], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &waves_[iw];
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
            component = std::complex<double>(components[is], components[is+1]);
            component /= sqrt(2.0);

            // Loop over waves in first star
            starPtr = &stars_[is];
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &waves_[iw];
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &stars_[is+1];
            component = conj(component);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &waves_[iw];
               if (!wavePtr->implicit) {
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
   void Basis<D>::convertFieldDftToComponents(RFieldDft<D>& dft, 
                                              DArray<double>& components)
   {
      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      Star* starPtr;                  // pointer to current star
      Wave* wavePtr;                  // pointer to current wave
      std::complex<double> component; // coefficient for star
      IntVec<D> indices;              // dft grid indices of wave
      int rank;                       // dft grid rank of wave
      int is;                         // star index

      // Loop over stars
      is = 0;
      while (is < nStar_) {
         starPtr = &stars_[is];
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Characteristic wave is first wave of star
            wavePtr = &waves_[starPtr->beginId];
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
            wavePtr = &waves_[starPtr->beginId];
            if (wavePtr->implicit) {
               starPtr = &stars_[is+1];
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &waves_[starPtr->endId-1];
               UTIL_CHECK(!(wavePtr->implicit));
            } 
            indices = wavePtr->indicesDft;
            rank = dftMesh.rank(indices);

            // Compute component value
            component = std::complex<double>(dft[rank][0], dft[rank][1]);
            UTIL_CHECK(abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);
            if (starPtr->invertFlag == -1) {
               component = conj(component);
            }
            components[is] = component.real();
            components[is+1] = component.imag();

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }

   template <int D>
   int Basis<D>::nBasis() const
   {
      int count = 0;
      for (int i = 0; i < stars_.capacity(); ++i) {
         if (!stars_[i].cancel) {
            count++;
         }
      }
      return count;
   }

   template <int D>
   bool Basis<D>::isValid() const
   {
      IntVec<D> v;
      int iw, is;

      // Loop over dft mesh to check consistency of waveIds_ and waves_
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         v = itr.position();
         iw = waveId(v);
         if (wave(iw).indicesDft != v) {
            std::cout << "Inconsistent waveId and Wave::indicesDft" 
                      << std::endl;
            return false;
         }
      }

      // Loop over waves, check consistency of wave data.
      for (iw = 0; iw < nWave_; ++iw) {

         // Check that indicesBz is an image of indicesDft
         v = waves_[iw].indicesBz;
         mesh().shift(v);
         if (v != waves_[iw].indicesDft) {
             std::cout << "shift(indicesBz) != indicesDft" << std::endl;
             return false;
         }
 
         // Check starId
         is = waves_[iw].starId;
         if (iw < stars_[is].beginId) {
             std::cout << "Value of Wave::starId < Star::beginId" << std::endl;
             return false;
         } 
         if (iw >= stars_[is].endId) {
             std::cout << "Value of Wave::starId >= Star::endId" << std::endl;
             return false;
         } 
      }

      // Loop over stars, check consistency of star data.
      for (is = 0; is < nStar_; ++is) {

         // Check star size
         if (stars_[is].size != stars_[is].endId - stars_[is].beginId) {
            std::cout << "Inconsistent Star::size" 
                      << std::endl;
             return false;
         }
         if (is > 0) {
            if (stars_[is].beginId != stars_[is-1].endId) {
               std::cout << "Star ranges not consecutive" << std::endl;
               return false;
            }
         }

         // Check star ids of waves in star
         for (iw = stars_[is].beginId; iw < stars_[is].endId; ++iw) {
            if (waves_[iw].starId != is) {
               std::cout << "Inconsistent Wave::starId and Star range" 
                         << std::endl;
               return false;
            }
         }

         // Check ordering of waves in star
         for (iw = stars_[is].beginId + 1; iw < stars_[is].endId; ++iw) {
            if (waves_[iw].indicesBz > waves_[iw-1].indicesBz) {
               std::cout << "Failure of ordering by indicesB within star" 
                         << std::endl;
               return false;
            }
            if (waves_[iw].indicesBz == waves_[iw-1].indicesBz) {
               std::cout << "Equal values of indicesBz within star" 
                         << std::endl;
               return false;
            }
         }

         // Loop over closed stars and related pairs of stars
         is = 0;
         while (is < nStar_) {
            if (stars_[is].invertFlag == 0) {
               ++is;
            } else {
               if (stars_[is].invertFlag != 1) {
                  std::cout << "Expected invertFlag == 1" << std::endl;
                  return false;
               }
               if (stars_[is+1].invertFlag != -1) {
                  std::cout << "Expected invertFlag == -1" << std::endl;
                  return false;
               }
               is += 2;
            }
         }

      }

      return true;
   }

}
}

#endif
