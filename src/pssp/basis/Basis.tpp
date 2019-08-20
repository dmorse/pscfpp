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
#include <fstream>

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
      nBasis_(0), 
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
         std::ifstream in;
         in.open(groupName);
         if (in.is_open()) {
           in >> group;
           UTIL_CHECK(group.isValid());
         } else {
           UTIL_THROW("Unknown space group");
         }
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

      // Apply validity test suite
      bool valid = isValid();
      if (!valid) {
         UTIL_THROW("Basis failed validity check");
      }

   }

   /*
   * Construct ordered list of waves.
   * 
   * On exit:
   *  - Array waves_ contains list of waves ordered by sqNorm.
   *  - Each wave has indicesDft, indicesBz and sqNorm set.
   *  - Array stars_  is still empty.
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
      * list - set of waves of equal norm, compared by indicesDft
      * star - set of symmetry-related waves, compared by indicesDft
      * tempStar - temporary star, sorted by descending indicesBz
      * tempList - temporary list, ordered by star
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
      double Gsq_max = 1.0E-8;
      double epsilon = 1.0E-8;
      double twoPi = 2.0*Constants::Pi;
      IntVec<D> rootVec;
      IntVec<D> vec;
      IntVec<D> nVec;
      int listId = 0;      // id for this list
      int listBegin = 0;   // id of first wave in this list
      int listEnd = 0;     // (id of last wave in this list) + 1
      int starId = 0;      // id for this star
      int starBegin = 0;   // id of first wave in this star
      int starSize = 0;    // size (number of waves) in this star
      int i, j, k;
      bool cancel;

      // Loop over all waves
      nBasis_ = 0;
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

         // Process completed list of wavectors of equal norm
         if (newList && listEnd > 0) {

            // Copy waves of equal norm into std::set "list"
            list.clear();
            tempList.clear();
            for (j = listBegin; j < listEnd; ++j) {
               wave.indicesDft = waves_[j].indicesDft;
               wave.indicesBz = waves_[j].indicesBz;
               wave.sqNorm = waves_[j].sqNorm;
               if (j > listBegin) {
                  UTIL_CHECK(abs(wave.sqNorm-waves_[j].sqNorm) < 2.0*epsilon); 
               }
               list.insert(wave);
            }

            // On entry to each iteration of the loop over stars,
            // rootIter and nextInvert are known. The iterator rootIter 
            // points to the wave that will be used as the root of the
            // next star. The flag nextInvert is equal to -1 iff the 
            // previous star was open under inversion, +1 otherwise.

            // Initial values for first star in this list
            rootItr = list.begin();
            int nextInvert = 1;

            // Loop over stars with a list of waves of equal norm
            while (list.size() > 0) {

               rootVec = rootItr->indicesBz;
               Gsq = unitCell().ksq(rootVec);
               cancel = false;

               // Construct a star from root vector, by applying every
               // symmetry operation in the group to the root wavevector.
               star.clear();
               for (j = 0; j < group.size(); ++j) {

                  // Apply symmetry (i.e., multiply by rotation matrix)
                  // vec = rotated wavevector.
                  vec = rootVec*group[j];

                  // Check that rotated vector has same norm as root.
                  UTIL_CHECK(abs(Gsq - unitCell().ksq(vec)) < epsilon);

                  // Initialize TWave object associated with rotated wave
                  wave.sqNorm = Gsq;
                  wave.indicesBz = vec;
                  wave.indicesDft = vec;
                  mesh().shift(wave.indicesDft);

                  // Compute phase for coeff. of wave in basis function.
                  // Convention -pi < phase <= pi.
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

                  // Check for cancellation of star: It can be shown
                  // that a star is cancelled iff application of any
                  // symmetry operation in the group yields a rotated
                  // vector vec equal to the root vector rootVec, but 
                  // with an nonzero phase, creating a contradiction.

                  if (vec == rootVec) {
                     if (abs(wave.phase) > 1.0E-6) {
                        cancel = true;
                     }
                  }

                  // Attempt to add wave to set star. Insertion fails 
                  // if star contains a TWave with the same indicesDft.
                  star.insert(wave);
               }
               starSize = star.size();

               // Append waves in star to the std::vector tempStar
               tempStar.clear();
               setItr = star.begin(); 
               for ( ; setItr != star.end(); ++setItr) {
                  tempStar.push_back(*setItr);
               }

               // Sort tempStar, ordered by indicesBz, in descending order
               // Duplicate values of indicesBz are permitted
               TWaveBzComp<D> waveBzComp;
               std::sort(tempStar.begin(), tempStar.end(), waveBzComp);
               
               // Append contents of tempStar to tempList, erase from list
               for (j = 0; j < tempStar.size(); ++j) {
                  list.erase(tempStar[j]);
                  tempList.append(tempStar[j]);
               }
               UTIL_CHECK(tempList.size()+list.size() == listEnd-listBegin);

               // If this star is not cancelled, increment the total
               // number of basis functions.
               if (!cancel) ++nBasis_;

               // Initialize a Star object (startInvert not yet known)
               newStar.eigen = Gsq;
               newStar.beginId = starBegin;
               newStar.endId = newStar.beginId + star.size();
               newStar.size = star.size();
               newStar.cancel = cancel;

               // Determine invertFlag, rootItr and nextInvert
               if (nextInvert == -1) {

                  // If this star is the 2nd of a pair related by symmetry,
                  // set root of next star to first wave of remaining list.

                  newStar.invertFlag = -1;
                  rootItr = list.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the second of a pair of partners,
                  // then determine if it is closed under inversion.

                  // Compute negation of root vector. 
                  nVec.negate(rootItr->indicesBz);

                  // Shift of the root to a DFT mesh
                  (*meshPtr_).shift(nVec);
       
                  // Search for negation of root vector within this star
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

                     // Negation should be found in star or remaining list
                     if (!negationFound) {
                        std::cout << "Negation not found for: " << "\n";
                        std::cout << " vec (ft):" << rootItr->indicesDft <<"\n"; 
                        std::cout << " vec (bz):" << rootItr->indicesBz <<"\n"; 
                        std::cout << "-vec (dft:" << nVec << "\n";
                     }
                     UTIL_CHECK(negationFound);
                  }

               }


               stars_.append(newStar);
               ++starId;
               starBegin = newStar.endId;
            } 
            UTIL_CHECK(list.size() == 0);
            UTIL_CHECK(tempList.size() == listEnd - listBegin);
            // Process of list is now complete

            // Copy tempList into corresponding section of waves_,
            // overwriting the section of waves_ used to create the list.
            // Compute a complex coefficient of unit norm for each wave.
            for (j = 0; j < tempList.size(); ++j) {
               k = j + listBegin;
               waves_[k].indicesDft = tempList[j].indicesDft;
               waves_[k].indicesBz = tempList[j].indicesBz;
               waves_[k].sqNorm = tempList[j].sqNorm;
               coeff = std::complex<double>(0.0, tempList[j].phase);
               waves_[k].coeff = exp(coeff);
            }
            // At this point, waves_[k].coeff has unit absolute magnitude,
            // and correct relative phases for waves within a star. The
            // root star has coefficient of unity.

            ++listId;
            listBegin = listEnd;
         }
      }
      nStar_ = stars_.size();

      std::complex<double> rootCoeff; // Coefficient of root wave
      std::complex<double> partCoeff; // Coefficient of partner of root
      std::complex<double> d;
      int rootId, partId;

      // Final processing of coefficients of waves in stars. Require
      // that the root wave of each star and its negation have conjugate 
      // coefficients, and each basis function is normalized.
      for (i = 0; i < nStar_; ++i) {

         // Treat open and closed stars separately

         if (stars_[i].invertFlag == 0) {
      
            rootId = stars_[i].beginId;
            stars_[i].waveBz = waves_[rootId].indicesBz;

            if (stars_[i].cancel) {

               // If the star is closed set all coefficients to zero
               std::complex<double> czero(0.0, 0.0);
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff = czero;
               }

            } else {

               // Compute negation of root vector, shift to DFT mesh
               nVec.negate(waves_[rootId].indicesBz);
               (*meshPtr_).shift(nVec);

               // Find negation of root in this star, set partId to index 
               bool negationFound = false;
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  if (nVec == waves_[j].indicesDft) {
                     partId = j;
                     negationFound = true;
                     break;
                  }
               }
               UTIL_CHECK(negationFound);

               // Divide all coefficients by the root coefficient
               rootCoeff = waves_[rootId].coeff;
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff /= rootCoeff;
               }

               if (partId !=  rootId) {

                  // Compute common divisor
                  partCoeff = waves_[partId].coeff;
                  d = sqrt(partCoeff);
                  if (abs(imag(d)) > 1.0E-8) {
                     if (imag(d) < 0.0) {
                        d = -d;
                     }
                  }

                  // Divide all coefficients by constant divisor
                  for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                     waves_[j].coeff /= d;
                  }

               }

            } // end if (cancel) ... else ...
   
         } else 
         if (stars_[i].invertFlag == 1) {

            // Process a pair of open stars related by inversion.
            // Preconditions:
            UTIL_CHECK(stars_[i].size == stars_[i+1].size);
            UTIL_CHECK(stars_[i].cancel == stars_[i+1].cancel);

            // Identify root of this star (star i)
            rootId = stars_[i].beginId;
            stars_[i].waveBz = waves_[rootId].indicesBz;

            // Compute negation of root vector, shift to DFT mesh
            nVec.negate(waves_[rootId].indicesBz);
            (*meshPtr_).shift(nVec);

            // Find negation of root wave in the next star (star i+1)
            bool negationFound = false;
            for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
               if (nVec == waves_[j].indicesDft) {
                  partId = j;
                  stars_[i+1].waveBz = waves_[j].indicesBz;
                  negationFound = true;
                  break;
               }
            }
            UTIL_CHECK(negationFound);

            if (stars_[i].cancel) {

               std::complex<double> czero(0.0, 0.0);
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff = czero;
               }
               for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
                  waves_[j].coeff = czero;
               }

            } else {

               // Divide all coefficients in this star by root coeff
               rootCoeff = waves_[rootId].coeff;
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff /= rootCoeff;
               }

               // Divide coefficients in next star by a partner coeff.
               partCoeff = waves_[partId].coeff;
               for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
                  waves_[j].coeff /= partCoeff;
               }
   
            } // end if (cancel) ... else ...

         } // end if (invertFlag == -1)
   
      } // end loop over stars

      // Final processing of waves in stars
      for (i = 0; i < nStar_; ++i) {

         // Set starId and normalize coefficients for associated waves
         double snorm = 1.0/sqrt(double(stars_[i].size));
         for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
            waves_[j].starId = i;
            waves_[j].coeff *= snorm;
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
            component = std::complex<double>(components[is], 
                                             -components[is+1]);
            component /= sqrt(2.0);

            // Loop over waves in first star
            starPtr = &stars_[is];
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &waves_[iw];
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &stars_[is+1];
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = conj(component);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &waves_[iw];
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
            components[is] = component.real();
            components[is+1] = -component.imag();

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }


   template <int D>
   int Basis<D>::nBasis() const
   {  return nBasis_; }


   template <int D>
   void Basis<D>::outputWaves(std::ostream& out) const
   {
      out << std::endl;
      out << "Waves:" << std::endl;
      int i, j;
      for (i = 0; i < nWave_; ++i) {
         out << Int(i,4);
         out << Int(waves_[i].starId, 4);
         out << " |";
         for (j = 0; j < D; ++j) {
            out << Int(waves_[i].indicesDft[j], 4);
         }
         out << " |";
         for (j = 0; j < D; ++j) {
            out << Int(waves_[i].indicesBz[j], 4);
         }
         out << " | ";
         out << Dbl(waves_[i].sqNorm, 12);
         out << "  " << Dbl(waves_[i].coeff.real(), 10);
         out << "  " << Dbl(waves_[i].coeff.imag(), 10);
         out << std::endl;
      }
   }


   template <int D>
   void Basis<D>::outputStars(std::ostream& out) const
   {
      out << std::endl;
      out << "Stars:" << std::endl;
      int i, j;
      for (i = 0; i < nStar_; ++i) {
         out << Int(i, 4)
             << Int(stars_[i].size, 3)
             << Int(stars_[i].beginId, 5)
             << Int(stars_[i].endId, 5)
             << Int(stars_[i].invertFlag, 3)
             << Int(stars_[i].cancel, 3);
         out << " |";
         for (j = 0; j < D; ++j) {
            out << Int(stars_[i].waveBz[j], 4);
         }
         out << " | " << Dbl(stars_[i].eigen, 12);
         out << std::endl;
      }
   }

   template <int D>
   bool Basis<D>::isValid() const
   {
      IntVec<D> v;
      int is, iw, iwp;

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

         // Check that wave indicesBz is an image of indicesDft
         v = waves_[iw].indicesBz;
         mesh().shift(v);
         if (v != waves_[iw].indicesDft) {
             std::cout << "shift(indicesBz) != indicesDft" << std::endl;
             return false;
         }
 
         // Check starId
         is = waves_[iw].starId;
         if (iw < stars_[is].beginId) {
             std::cout << "Wave::starId < Star::beginId" << std::endl;
             return false;
         } 
         if (iw >= stars_[is].endId) {
             std::cout << " Wave::starId >= Star::endId" << std::endl;
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
         std::complex<double> cdel;
         int begin, end;
         bool negationFound, cancel;
         is = 0;
         while (is < nStar_) {
            cancel = stars_[is].cancel;

            if (stars_[is].invertFlag == 0) {
            
               // Test star that is closed under inversion

               begin = stars_[is].beginId; 
               end = stars_[is].endId; 
               for (iw = begin; iw < end; ++iw) {
                  v.negate(waves_[iw].indicesBz);
                  (*meshPtr_).shift(v);
                  negationFound = false;
                  for (iwp = begin; iw < end; ++iwp) {
                     if (waves_[iwp].indicesDft == v) {
                        negationFound = true;
                        if (!cancel) {
                           cdel = conj(waves_[iwp].coeff);
                           cdel -= waves_[iw].coeff;
                        }
                        break;
                     }
                  }
                  if (!negationFound) {
                     std::cout << "Negation not found in closed star" 
                               << std::endl;
                     return false;
                  }
                  if (!cancel && abs(cdel) > 1.0E-8) {
                     std::cout << "Coefficients not conjugates in closed star" 
                               << std::endl;
                     return false;
                  }
               }
 
               ++is;

            } else {

               // Test pairs of open stars

               if (stars_[is].invertFlag != 1) {
                  std::cout << "Expected invertFlag == 1" << std::endl;
                  return false;
               }
               if (stars_[is+1].invertFlag != -1) {
                  std::cout << "Expected invertFlag == -1" << std::endl;
                  return false;
               }
               if (stars_[is+1].size != stars_[is].size) {
                  std::cout << "Parners of different size" << std::endl;
                  return false;
               }

               begin = stars_[is+1].beginId; 
               end = stars_[is+1].endId;

               // Check existence of negation and conjugate coefficients
               // Loop over waves in first star
               for (iw = stars_[is].beginId; iw < stars_[is].endId; ++iw) {
                  v.negate(waves_[iw].indicesBz);
                  (*meshPtr_).shift(v);
                  negationFound = false;
                  // Loop over second star, searching for negation
                  for (iwp = begin; iw < end; ++iwp) {
                     if (waves_[iwp].indicesDft == v) {
                        negationFound = true;
                        if (!cancel) {
                           cdel = conj(waves_[iwp].coeff);
                           cdel -= waves_[iw].coeff;
                           if (abs(cdel) > 1.0E-8) {
                              std::cout  <<
                                 "Coefficients not conjugates for open star" 
                                 << std::endl;
                              return false;
                           }
                        }
                        break;
                     }
                  }
                  if (!negationFound) {
                     std::cout << "Negation not found for open star" 
                               << std::endl;
                     return false;
                  }
               }

               is += 2;
            }
         }

      }
 
      // If end of function is reached, then all tests passed
      return true;
   }

}
}

#endif
