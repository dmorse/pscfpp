#ifndef PSSP_BASIS_TPP
#define PSSP_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"
#include "TWave.h"
#include "groupFile.h"
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <algorithm>
#include <vector>
#include <set>
#include <fstream>

namespace Pscf {

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
   * Destructor.
   */
   template <int D>
   Basis<D>::~Basis()
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
         // Create identity group by default
         group.makeCompleteGroup();
      } else {
         bool foundFile = false;
         {
            std::ifstream in;
            in.open(groupName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
               foundFile = true;
            }
         }
         if (!foundFile) {
            std::string fileName = makeGroupFileName(D, groupName);
            std::ifstream in;
            in.open(fileName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
            } else {
               Log::file() << "\nFailed to open group file: " 
                           << fileName << "\n";
               Log::file() << "\n Error: Unknown space group\n";
               UTIL_THROW("Unknown space group");
            }
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
      double Gsq_max;
      double phase_diff;
      const double twoPi = 2.0*Constants::Pi;
      const double epsilon = 1.0E-8;
      IntVec<D> meshDimensions = mesh().dimensions();
      IntVec<D> rootVecBz;   // BZ indices for root of this star
      IntVec<D> vec;         // Indices of temporary wavevector
      IntVec<D> nVec;        // Indices of negation of a wavevector
      int listId = 0;        // id for this list
      int listBegin = 0;     // id of first wave in this list
      int listEnd = 0;       // (id of last wave in this list) + 1
      int listSize;          // listEnd - listBegin
      int starId = 0;        // id for this star
      int starBegin = 0;     // id of first wave in this star
      int i, j, k;
      bool cancel;

      /*
      * Overview of algorithm:
      *
      * Loop over index i of array waves_ {
      *
      *   Search for end of a "list", contiguous block of waves of equal
      *   magnitude, with indices [listBegin,listEnd-1]. Set newList true.
      *
      *   // Each list can contain one or more stars.
      *   // Process the newly identified list to identify stars
      *
      *   If (newList) {
      *   
      *     Copy all waves in range into container std::set<TWave<D>> list
      *
      *     Set rootItr to first wave in list
      *
      *     // Loop over stars within list
      *     while (list.size() > 0) {
      *
      *       // To generate a star from a root wave rootItr, 
      *       // loop over symmetry operations of space group.
      *       For each group symmetry operation group[j] {
      *         Compute vec = (rootItr->indicesBz)*group[j]
      *         Set phase = rootItr->indicesBz .dot. group[j].t 
      *         Check for cancellation of the star, set cancel flag
      *         Add wave to std::set<TWave> star if not added before
      *         // Here, use of a std::set simplifies test of uniqueness
      *       }
      *
      *       Copy all waves from star to std::vector<TWave> tempStar
      *       Sort tempStar by indicesBz in descending order
      *       // Here, use of ordered container std::vector allows sorting
      *       For each wave in tempStar {
      *         Append wave to std::vector<TWave> tempList
      *         Erase wave from std::set<TWave> list
      *         // Here, use of a std::set for list simplifies erasure
      *       }
      *
      *       Initialize a Star object newStar and assign values
      *       to members beginId, endId, size, eigen, cancel
      *
      *       // Assign values of newStar.invertFlag, rootItr, nextInvert
      *       if (nextInvert == -1) { 
      *          // This is the second star in pair
      *          newStar.invertFlag = -1;
      *          nextInvert = 1;
      *          Set rootItr to the first wave in remaining list
      *       } else {
      *          Search for negation of rootItr in this star
      *          if negation is in this star {
      *             newStar.invertFlag = 0
      *             nextInvert = 1;
      *             Set rootItr to the first wave in remaining list
      *          } else 
      *          Search for negation of rootItr in this remaining list
      *          if the negation is in the remaining list {
      *             newStar.invertFlag = 1
      *             nextInvert = -1;
      *             set rootItr to negation of current root
      *          }
      *       }
      *
      *       Append newStar object to DArray<Basis::Star> stars_ 
      *
      *     } // end loop over stars in a single list
      *
      *     For each TWave in tempList {
      *       Copy TWave in tempList to a Basis:Wave in  waves_
      *       Assign a complex coefficient of unit norm to the Wave
      *     }
      *     // This overwrites a block of the waves_ array 
      *
      *     // At this point, coefficients of waves have correct
      *     // correct relative phases within a star, but not final 
      *     // absolute phases and have unit absolute magnitude.
      *
      *   } // end processing of one list (waves of equal norm)
      *
      * } // End initial processing of all waves and stars
      *
      * // Set absolute wave coefficients
      * For each star in array stars_ {
      *   if star is closed (star.invertFlag == 0) {
      *     if star is cancelled {
      *       set all coefficients to zero
      *     } else {
      *       Set the root to the first wave in the star
      *       Check that negation of root is also in this star
      *       Divide all coefficients by the root coefficient
      *       Require coeffs of root & negation be complex conjugates
      *     }
      *   } else
      *   if (star.invertFlag == 1) {
      *     Set root of this star to the 1st wave in the star
      *     Find negation of root (aka "partner") in next star.
      *     If star is cancelled {
      *       Set coefficients in this star and next to zero
      *     } else {
      *       Divide coefficients in this star by root coefficient
      *       Divide coefficients in next star by partner coefficient
      *     }
      *   }
      *   // Note: If star.invertFlag = -1, do nothing
      * }
      *
      * // For all waves, normalize coefficients and set starId
      * For each star in stars_ {
      *   For each wave in this star {
      *     Set Wave.starId
      *     Divide coefficient by sqrt(double(star.size))
      *   }
      * }
      *
      * // For all waves, set implicit member and add to look up table
      * For each wave in waves_ { 
      *   Set implicit attribute
      *   Assign waveIds_[rank] = i
      * }
      */

      // Loop over all waves
      nBasis_ = 0;
      nBasisWave_ = 0;
      Gsq_max = waves_[0].sqNorm;
      for (i = 1; i <= nWave_; ++i) {

         // Determine if this wave begins a new list
         bool newList = false;
         if (i == nWave_) {
            listEnd = i;
            listSize = listEnd - listBegin;
            newList = true;
         } else {
            Gsq = waves_[i].sqNorm;
            if (Gsq > Gsq_max + epsilon) {
               Gsq_max = Gsq;
               listEnd = i;
               listSize = listEnd - listBegin;
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
                  UTIL_CHECK( abs(wave.sqNorm-waves_[j].sqNorm) 
                                 < 2.0*epsilon ); 
               }
               list.insert(wave);
            }

            // On entry to each iteration of the loop over stars,
            // rootItr and nextInvert are known. The iterator rootItr 
            // points to the wave in the remaining list that will be 
            // used as the root of the next star. The flag nextInvert 
            // is equal to -1 iff the previous star was the first of
            // a pair that are open under inversion, and is equal 
            // to + 1 otherwise.

            // Initial values for first star in this list
            rootItr = list.begin();
            int nextInvert = 1;

            // Loop over stars with a list of waves of equal norm,
            // removing each star from the list as it is identified.
            // The root of the next star must have been chosen on
            // entry to each iteration of this loop.

            while (list.size() > 0) {

               rootVecBz = rootItr->indicesBz;
               Gsq = rootItr->sqNorm;
               cancel = false;
               star.clear();

               // Construct a star from root vector, by applying every
               // symmetry operation in the group to the root wavevector.
               for (j = 0; j < group.size(); ++j) {

                  // Apply symmetry (i.e., multiply by rotation matrix)
                  // vec = rotated wavevector.
                  vec = rootVecBz*group[j];

                  // Check that rotated vector has same norm as root.
                  UTIL_CHECK(abs(Gsq - unitCell().ksq(vec)) < epsilon);

                  // Initialize TWave object associated with rotated wave
                  wave.sqNorm = Gsq;
                  wave.indicesBz = shiftToMinimum(vec, meshDimensions, 
                                                  *unitCellPtr_);
                  wave.indicesDft = vec;
                  mesh().shift(wave.indicesDft);

                  // Compute phase for coeff. of wave in basis function.
                  // Convention -pi < phase <= pi.
                  wave.phase = 0.0;
                  for (k = 0; k < D; ++k) {
                     wave.phase += rootVecBz[k]*(group[j].t(k));
                  }
                  while (wave.phase > 0.5) {
                     wave.phase -= 1.0;
                  }
                  while (wave.phase <= -0.5) {
                     wave.phase += 1.0;
                  }
                  wave.phase *= twoPi;

                  // Check for cancellation of star: The star is
                  // cancelled if application of any symmetry operation
                  // in the group to the root vector yields a rotated 
                  // vector equivalent to the root vector but with a 
                  // nonzero phase, creating a contradiction.

                  if (wave.indicesDft == rootItr->indicesDft) {
                     if (abs(wave.phase) > 1.0E-6) {
                        cancel = true;
                     }
                  }

                  // Search for an equivalent wave already in the star.
                  // Note: Equivalent waves have equal DFT indices.
                  setItr = star.find(wave);

                  if (setItr == star.end()) {

                     // If no equivalent wave is found in the star,
                     // then add this wave to the star
                     star.insert(wave);

                  } else {

                     // If an equivalent wave is found, check if the
                     // phases are equivalent. If not, the star is
                     // cancelled.

                     phase_diff = setItr->phase - wave.phase;
                     while (phase_diff > 0.5) {
                        phase_diff -= 1.0;
                     }
                     while (phase_diff <= -0.5) {
                        phase_diff += 1.0;
                     }
                     if (abs(phase_diff) > 1.0E-6) {
                        cancel = true;
                     }

                  }

               }

               // Copy all waves from set star to std::vector tempStar
               tempStar.clear();
               setItr = star.begin(); 
               for ( ; setItr != star.end(); ++setItr) {
                  tempStar.push_back(*setItr);
               }

               // Sort tempStar, in descending order by indicesBz.
               TWaveBzComp<D> waveBzComp;
               std::sort(tempStar.begin(), tempStar.end(), waveBzComp);
               
               // Append contents of tempStar to tempList, erase from list
               int tempStarSize = tempStar.size();
               for (j = 0; j < tempStarSize; ++j) {
                  list.erase(tempStar[j]);
                  tempList.append(tempStar[j]);
               }
               UTIL_CHECK((int)(tempList.size()+list.size()) == listSize);

               // If this star is not cancelled, increment the number of 
               // basis functions (nBasis_) and waves in basis (nBasisWave_)
               if (!cancel) {
                  ++nBasis_;
                  nBasisWave_ += star.size();
               }

               // Initialize a Star object 
               newStar.eigen = Gsq;
               newStar.beginId = starBegin;
               newStar.endId = newStar.beginId + star.size();
               newStar.size = star.size();
               newStar.cancel = cancel;
               // Note: newStar.starInvert is not yet known

               // Determine invertFlag, rootItr and nextInvert
               if (nextInvert == -1) {

                  // If this star is 2nd of a pair related by inversion,
                  // set root of next star to 1st wave of remaining list.

                  newStar.invertFlag = -1;
                  rootItr = list.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the 2nd of a pair of partners,
                  // then determine if it is closed under inversion.

                  // Compute negation nVec of root vector in FBZ
                  nVec.negate(rootVecBz);

                  // Shift negation nVec to a DFT mesh
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

                     // If this star is closed under inversion, the root
                     // of next star is the 1st vector of remaining list.

                     newStar.invertFlag = 0;
                     rootItr = list.begin();
                     nextInvert = 1;

                  } else {

                     newStar.invertFlag = 1;
                     nextInvert = -1;

                     // If star is not closed, find negation of root in
                     // the remaining list, and use this negation as the
                     // root of the next star.

                     setItr = list.begin();
                     for ( ; setItr != list.end(); ++setItr) {
                        if (nVec == setItr->indicesDft) {
                           negationFound = true;
                           rootItr = setItr;
                           break;
                        }
                     }
                     // If negationFound, then rootItr->indicesDft = nVec

                     // Failure to find the negation here is an error:
                     // It must be either in this star or remaining list

                     if (!negationFound) {
                        std::cout << "Negation not found for: " << "\n";
                        std::cout << " vec (ft):" 
                                  << rootItr->indicesDft <<"\n"; 
                        std::cout << " vec (bz):" 
                                  << rootItr->indicesBz <<"\n"; 
                        std::cout << "-vec (dft):" << nVec << "\n";
                        UTIL_CHECK(negationFound);
                     }

                  }

               }

               stars_.append(newStar);
               ++starId;
               starBegin = newStar.endId;

            } 
            // End loop over stars within a list.

            UTIL_CHECK(list.size() == 0);
            UTIL_CHECK(tempList.size() == listEnd - listBegin);

            // Copy tempList into corresponding section of waves_,
            // overwriting the section of waves_ used to create the list.
            // Compute a complex coefficient of unit norm for each wave.
            for (j = 0; j < tempList.size(); ++j) {
               k = j + listBegin;
               waves_[k].indicesDft = tempList[j].indicesDft;
               waves_[k].indicesBz = tempList[j].indicesBz;
               waves_[k].sqNorm = tempList[j].sqNorm;
               coeff = std::complex<double>(0.0, tempList[j].phase);
               coeff = exp(coeff);
               if (abs(imag(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(real(coeff), 0.0);
               }
               if (abs(real(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(0.0, imag(coeff));
               }
               waves_[k].coeff = coeff;
            }

            // Processing of list is now complete.
            // Here, waves_[k].coeff has unit absolute magnitude, and
            // correct relative phases for waves within a star, but
            // the coeff may not be unity for the first or last wave in
            // the star.

            ++listId;
            listBegin = listEnd;
         } 
         // Finished processing a list of waves of equal norm

      } // End loop over all waves
      nStar_ = stars_.size();
      // Complete initial processing of all lists and stars

      std::complex<double> rootCoeff; // Coefficient of root wave
      std::complex<double> partCoeff; // Coefficient of partner of root
      std::complex<double> d;
      int rootId, partId;

      // Final processing of coefficients of waves in stars. 
      // Require that the root wave of each star and its negation 
      // have conjugate coefficients, and that each basis function 
      // is normalized.
      for (i = 0; i < nStar_; ++i) {

         // Treat open and closed stars separately

         if (stars_[i].invertFlag == 0) {
     
            // Identify root of this star (star i)
            // Set the root to be the first wave in the star

            rootId = stars_[i].beginId;
            stars_[i].waveBz = waves_[rootId].indicesBz;

            if (stars_[i].cancel) {

               // If the star is cancelled, set all coefficients to zero
               std::complex<double> czero(0.0, 0.0);
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff = czero;
               }

            } else { // if not cancelled

               // Compute nVec = negation of root, shifted to DFT mesh
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

               // For invertFlag == 0, failure to find nVec in this
               // star is a fatal error
               UTIL_CHECK(negationFound);

               // Divide all coefficients by the root coefficient
               rootCoeff = waves_[rootId].coeff;
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff /= rootCoeff;
               }

               // Require coefficients of root and negation are conjugates
               if (partId != rootId) {

                  // Compute common divisor
                  partCoeff = waves_[partId].coeff;
                  d = sqrt(partCoeff);
                  if (abs(imag(d)) > 1.0E-8) {
                     if (imag(d) < 0.0) {
                        d = -d;
                     }
                     for (j=stars_[i].beginId; j < stars_[i].endId; ++j) {
                        waves_[j].coeff /= d;
                     }
                  }

               }

            } // end if (cancel) ... else ...
   
            // end if (stars_[i].invertFlag == 0) 
         } else 
         if (stars_[i].invertFlag == 1) {

            // Process a pair of open stars related by inversion.

            // Preconditions:
            UTIL_CHECK(stars_[i].size == stars_[i+1].size);
            UTIL_CHECK(stars_[i].cancel == stars_[i+1].cancel);
            // UTIL_CHECK(stars_[i+1].invertFlag == -1);

            // Identify root of this star (star i)
            // Set the root to be the first wave in the star
            rootId = stars_[i].beginId;
            stars_[i].waveBz = waves_[rootId].indicesBz;

            // Compute nVec = negation of root vector, shifted to DFT mesh
            nVec.negate(waves_[rootId].indicesBz);
            (*meshPtr_).shift(nVec);

            // Seek negation of root wave in the next star (star i+1)
            bool negationFound = false;
            for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
               if (nVec == waves_[j].indicesDft) {
                  partId = j;
                  stars_[i+1].waveBz = waves_[j].indicesBz;
                  negationFound = true;
                  break;
               }
            }

            // For invertFlag == 1, absence of nVec in next star is fatal
            UTIL_CHECK(negationFound);

            if (stars_[i].cancel) {

               std::complex<double> czero(0.0, 0.0);
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff = czero;
               }
               for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
                  waves_[j].coeff = czero;
               }

            } else { // if star is not cancelled

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

         }  // end if (invertFlag==0) ... else if (invertFlag==1) ...

         // Note: If invertFlag == -1, do nothing and continue
         // Related stars with invertFlag == 1 and -1 are treated together
   
      } // end loop over stars

      // For all waves, normalize coefficients and set starId
      for (i = 0; i < nStar_; ++i) {
         double snorm = 1.0/sqrt(double(stars_[i].size));
         for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
            waves_[j].coeff *= snorm;
            waves_[j].starId = i;
         }
      }

      // For each wave, set implicit attribute and add to look-up table
      for (i = 0; i < nWave_; ++i) {
         vec = waves_[i].indicesDft;

         // Validity check - check that vec is in dft mesh
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

 
   /*
   *  Update wave norms after change in unit cell dimensions.
   */
   template <int D>
   void Basis<D>::update()
   {
      IntVec<D> vec;

      // Process stars
      for (int i = 0; i < nStar_; ++i) {
         vec = stars_[i].waveBz;
         stars_[i].eigen = unitCell().ksq(vec);
      }

      // Process waves
      for (int i = 0; i < nWave_; ++i) {
         vec = waves_[i].indicesBz;
         waves_[i].sqNorm = unitCell().ksq(vec);
      }

   }

   // Return value of nBasis
   template <int D>
   int Basis<D>::nBasis() const
   {  return nBasis_; }


   template <int D>
   void Basis<D>::outputWaves(std::ostream& out, bool outputAll) const
   {
      out << "N_wave" << std::endl;
      if (outputAll) {
         out << "             " << nWave_ << std::endl;
      } else {
         out << "             " << nBasisWave_ << std::endl;
      }
      int i, j, k, starId;
      k = 0;
      for (i = 0; i < nWave_; ++i) {
         starId = waves_[i].starId;
         if (outputAll || (!stars_[starId].cancel)) {
            out << Int(k, 8);
            // out << " |";
            // for (j = 0; j < D; ++j) {
            //   out << Int(waves_[i].indicesDft[j], 4);
            // }
            // out << " |";
            for (j = 0; j < D; ++j) {
               out << Int(waves_[i].indicesBz[j], 5);
            }
            // out << " | ";
            out << Int(waves_[i].starId, 6);
            out << "  " << Dbl(waves_[i].coeff.real(), 15);
            out << "  " << Dbl(waves_[i].coeff.imag(), 15);
            out << Dbl(waves_[i].sqNorm, 20);
            out << std::endl;
            k++;
         }
      }
   }


   template <int D>
   void Basis<D>::outputStars(std::ostream& out, bool outputAll) const
   {
      // Output number of stars in appropriate format
      if (outputAll) {
          out << "N_star" << std::endl
              << "                 " << nStar_ << std::endl;
      } else {
          out << "N_basis" << std::endl
              << "                 " << nBasis_ << std::endl;
      }

      // Loop over stars
      int i, j, k;
      k = 0;
      for (i = 0; i < nStar_; ++i) {
         if (outputAll || (!stars_[i].cancel)) {
            out << Int(k, 5)
                << Int(stars_[i].size, 5)
                << Int(stars_[i].beginId, 8)
                << Int(stars_[i].endId, 8)
                << Int(stars_[i].invertFlag, 4);
            if (outputAll) {
               out << Int(stars_[i].cancel, 4);
            }
            //out << " |";
            for (j = 0; j < D; ++j) {
               out << Int(stars_[i].waveBz[j], 6);
            }
            //out << " | " 
            out << Dbl(stars_[i].eigen, 15);
            out << std::endl;
            ++k;
         }
      }
   }

   template <int D>
   bool Basis<D>::isValid() const
   {
      double Gsq;
      IntVec<D> v;
      int is, iw, iwp, j;

      // Check total number of waves == # of grid points
      if (nWave_ != mesh().size()) {
         std::cout << "nWave != size of mesh" << std::endl;
      }

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

      // Loop over elements of waves_, check consistency of wave data.
      for (iw = 0; iw < nWave_; ++iw) {

         // Check sqNorm
         v = waves_[iw].indicesBz;
         Gsq = unitCell().ksq(v);
         if (abs(Gsq - waves_[iw].sqNorm) > 1.0E-8) {
             std::cout << "Incorrect sqNorm:" << "\n"
                       << "wave.indicesBz = " << "\n"
                       << "wave.sqNorm    = " << waves_[iw].sqNorm << "\n"
                       << "|v|^{2}        = " << Gsq << "\n";
         }
 
         // Check that wave indicesBz is an image of indicesDft
         mesh().shift(v);
         if (v != waves_[iw].indicesDft) {
             std::cout << "shift(indicesBz) != indicesDft" << std::endl;
             return false;
         }

         // Compare Wave::starId to Star::beginId and Star::endId
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

      // Loop over all stars (elements of stars_ array)
      int nWave = 0;
      for (is = 0; is < nStar_; ++is) {

         // Check star size
         nWave += stars_[is].size;
         if (stars_[is].size != stars_[is].endId - stars_[is].beginId) {
            std::cout << "Inconsistent Star::size:" << std::endl;
            std::cout << "Star id    "  << is << std::endl;
            std::cout << "star size  "  << stars_[is].size << std::endl;
            std::cout << "Star begin "  << stars_[is].beginId << std::endl;
            std::cout << "Star end   "  << stars_[is].endId << std::endl;
            return false;
         }
         if (is > 0) {
            if (stars_[is].beginId != stars_[is-1].endId) {
               std::cout << "Star ranges not consecutive:" << std::endl;
               std::cout << "Star id    "     << is << std::endl;
               std::cout << "stars_[" << is << "]"   << ".beginId = " 
                         << stars_[is].beginId << std::endl;
               std::cout << "stars_[" << is-1 << "]" << ".endId   = " 
                         << stars_[is-1].endId << std::endl;
               return false;
            }
         }

         // Check star ids of waves in star
         for (iw = stars_[is].beginId; iw < stars_[is].endId; ++iw) {
            if (waves_[iw].starId != is) {
               std::cout << "Inconsistent Wave::starId :" << std::endl;
               std::cout << "star id      "  << is << std::endl;
               std::cout << "star beginId "  << stars_[is].beginId << "\n";
               std::cout << "star endId   "  << stars_[is].endId << "\n";
               std::cout << "wave id      "  << iw << "\n";
               std::cout << "wave starId  "  << waves_[iw].starId << "\n";
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

      } // End do loop over all stars

      // Check that all waves in mesh are accounted for in stars
      if (stars_[nStar_-1].endId != mesh().size()) {
         std::cout << "Star endI of last star != mesh size" << std::endl;
      }
      if (nWave != mesh().size()) {
         std::cout << "Sum of star sizes != mesh size" << std::endl;
      }

      // Loop over closed stars and related pairs of stars.
      // Test closure under inversion and conjugacy of coefficients.
      std::complex<double> cdel;
      bool negationFound, cancel;
      is = 0;
      while (is < nStar_) {
         cancel = stars_[is].cancel;

         if (stars_[is].invertFlag == 0) {
         
            // Test that star is closed under inversion and real
            int begin = stars_[is].beginId;
            int end = stars_[is].endId;
            for (iw = begin; iw < end; ++iw) {
               v.negate(waves_[iw].indicesBz);
               mesh().shift(v);
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
                  std::cout << "G = " << waves_[iw].indicesBz
                            << "coeff = " << waves_[iw].coeff 
                            << std::endl;
                  std::cout << "All waves in star " << is << "\n";
                  for (j=begin; j < end; ++j) {
                     std::cout << waves_[j].indicesBz << "  "
                               << waves_[j].coeff << "\n";
                  }
                  return false;
               }
               if (!cancel && abs(cdel) > 1.0E-8) {
                  std::cout << "Function for closed star is not real:" 
                            << "\n";
                  std::cout << "+G = " << waves_[iw].indicesBz
                            << "  coeff = " << waves_[iw].coeff 
                            << "\n";
                  std::cout << "-G = " << waves_[iwp].indicesBz
                            << "  coeff = " << waves_[iwp].coeff 
                            << "\n";
                  std::cout << "Coefficients are not conjugates." << "\n";
                  std::cout << "All waves in star " << is << "\n";
                  for (j=begin; j < end; ++j) {
                     std::cout << waves_[j].indicesBz << "  "
                               << waves_[j].coeff << "\n";
                  }
                  return false;
               }
               if (cancel && abs(waves_[iw].coeff) > 1.0E-8) {
                  std::cout << "Nonzero coefficient in a cancelled star" 
                            << "\n";
                  std::cout << "G = " << waves_[iw].indicesBz
                            << "  coeff = " << waves_[iw].coeff 
                            << "\n";
                  return false;
               }
            }

            // Finished processing a closed star, increment counter is
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
               std::cout << "Partner stars of different size" << std::endl;
               return false;
            }
            if (stars_[is+1].cancel != stars_[is].cancel) {
               std::cout << "Partners stars with different cancel flags" 
                         << std::endl;
               return false;
            }

            // Begin and end wave ids for the first and second stars
            int begin1 = stars_[is].beginId; 
            int end1 = stars_[is].endId;
            int begin2 = stars_[is+1].beginId; 
            int end2 = stars_[is+1].endId;

            // Check existence of negation and conjugate coefficients
            // Loop over waves in first star
            for (iw = begin1; iw < end1; ++iw) {
               v.negate(waves_[iw].indicesBz);
               mesh().shift(v);
               negationFound = false;
               // Loop over second star, searching for negation
               for (iwp = begin2; iw < end2; ++iwp) {
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
                  std::cout << "Negation not found for G in open star" 
                            << std::endl;
                  std::cout << "First star id = " << is << std::endl;
                  std::cout << "+G = " << waves_[iw].indicesBz
                            << "coeff = " << waves_[iw].coeff 
                            << std::endl;
                  std::cout << "Waves in star " << is 
                            << "  (starInvert ==1):" << "\n";
                  for (j = begin1; j < end1; ++j) {
                     std::cout << waves_[j].indicesBz  << "  "
                               << waves_[j].coeff << "\n";
                  }
                  std::cout << "Waves in star " << is+1 
                            << "  (starInvert == -1):" << "\n";
                  for (j=begin2; j < end2; ++j) {
                     std::cout << waves_[j].indicesBz  << "  "
                               << waves_[j].coeff << "\n";
                  }
                  return false;
               } else 
               if (!cancel && abs(cdel) > 1.0E-8) {
                  std::cout << "Error of coefficients in open stars:" 
                            << "\n";
                  std::cout << "First star id = " << is << std::endl;
                  std::cout << "+G = " << waves_[iw].indicesBz
                            << "  coeff = " << waves_[iw].coeff 
                            << "\n";
                  std::cout << "-G = " << waves_[iwp].indicesBz
                            << "  coeff = " << waves_[iwp].coeff 
                            << "\n";
                  std::cout << "Coefficients are not conjugates." 
                            << "\n";
                  std::cout << "Waves in star " << is 
                            << "  (starInvert ==1):" << "\n";
                  for (j = begin1; j < end1; ++j) {
                     std::cout << waves_[j].indicesBz  << "  "
                               << waves_[j].coeff << "\n";
                  }
                  std::cout << "Waves in star " << is+1 
                            << "  (starInvert == -1):" << "\n";
                  for (j=begin2; j < end2; ++j) {
                     std::cout << waves_[j].indicesBz  << "  "
                               << waves_[j].coeff << "\n";
                  }
                  return false;
               }
            }

            // Finished processing a pair, increment star counter by 2
            is += 2;

         } // end if (stars_[is].invertFlag == 0) ... else ...

      } // end while (is < nStar_) loop over stars
 
      // The end of this function is reached iff all tests passed.
      return true;
   }

}
#endif
