#ifndef PRDC_BASIS_TPP
#define PRDC_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"
#include "TWave.h"
#include "groupFile.h"
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/signal/Signal.h>

#include <algorithm>
#include <vector>
#include <set>
#include <fstream>

namespace Pscf {
namespace Prdc {

   /*
   * Constructor.
   */
   template <int D>
   Basis<D>::Basis()
    : waves_(),
      stars_(),
      waveIds_(),
      starIds_(),
      nWave_(0),
      nBasisWave_(0),
      nStar_(0),
      nBasis_(0),
      signalPtr_(nullptr),
      unitCellPtr_(0),
      meshPtr_(0),
      isInitialized_(false)
   {
      signalPtr_ = new Signal<void>();
   }

   /*
   * Destructor.
   */
   template <int D>
   Basis<D>::~Basis()
   {
      delete signalPtr_;
   }

   /*
   * Construct basis for pseudo-spectral scft.
   */
   template <int D>
   void Basis<D>::makeBasis(Mesh<D> const & mesh,
                            UnitCell<D> const & unitCell,
                            std::string groupName)
   {
      SpaceGroup<D> group;
      readGroup(groupName, group);
      makeBasis(mesh, unitCell, group);
   }

   /*
   * Construct a symmetry-adapted basis for pseudo-spectral scft.
   */
   template <int D>
   void Basis<D>::makeBasis(Mesh<D> const & mesh,
                            UnitCell<D> const & unitCell,
                            SpaceGroup<D> const & group)
   {
      // Precondition: Check compatibility of mesh with space group
      group.checkMeshDimensions(mesh.dimensions());

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
         UTIL_THROW("Basis failed validity test suite");
      }

      // Mark as initialized
      isInitialized_ = true;

      // Notify any observers of successful basis initialization
      signal().notify();
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

   /*
   * Complete construction of a basis by grouping presorted waves into
   * stars and completing initialization of all Wave and Star objects.
   */
   template <int D>
   void Basis<D>::makeStars(SpaceGroup<D> const & group)
   {

      /*
      * Conceptual definitions:
      *
      *   - A "list" is a set of wavevectors of equal magnitude.
      *   - A "star" is a set of wavevectors that are related by symmetry.
      *
      * Each list may contain one or more complete stars. Complete lists 
      * are identified as an intermediate step in identification of stars.
      *
      * During initial processing, wavevectors are temporarily stored in 
      * TWave<D> objects.  The following local containers of TWave<D> 
      * objects are used:
      *
      *   list - a std::set of waves of equal norm (a "list")
      *   star - a std::set of symmetry-related waves (a "star")
      *   tempStar - a sorted star, sorted by descending indicesBz
      *   tempList - a sorted list, with contiguous sorted stars
      *
      * Reasons for choice of some C++ standard lib container types:
      *  std::set for "list", "star" prevents duplicates, allows removal
      *  std::vector for tempStar allows use of std::sort 
      */

      // Local TWave<D> containers and associated iterators
      std::set< TWave<D>, TWaveDftComp<D> > list;
      std::set< TWave<D>, TWaveDftComp<D> > star;
      std::vector< TWave<D> > tempStar;
      GArray< TWave<D> > tempList;
      typename std::set< TWave<D>, TWaveDftComp<D> >::iterator rootItr;
      typename std::set< TWave<D>, TWaveDftComp<D> >::iterator setItr;

      // Local variables
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
      IntVec<D> rootVecDft;  // DFT indices for root of this star
      IntVec<D> vec;         // Indices of temporary wavevector
      IntVec<D> nVec;        // Indices of inverse of a wavevector
      int listBegin = 0;     // id of first wave in this list
      int listEnd = 0;       // (id of last wave in this list) + 1
      int listSize;          // listEnd - listBegin
      int starBegin = 0;     // id of first wave in this star
      int i, j, k;
      bool cancel;

      /*
      * Overview of algorithm:
      *
      * Precondition: Wavevectors in the array waves_ are sorted in
      * nondecreasing order by wavevector norm.
      *
      * Loop over index i of array waves_ {
      *
      *   Search for end of a "list" (i.e., contiguous block of waves
      *   of equal magnitude) by identifying changes in magnitude.
      *   The resulting list has indices [listBegin,listEnd-1].
      *   Set newList true.
      *
      *   // Each list may contain one or more stars.
      *   // Process the newly identified list to identify stars
      *   If (newList) {
      *
      *     Copy all waves in the range into std::set list
      *
      *     Set rootItr to the first wave in the list
      *
      *     // Loop over stars within the list
      *     while (list.size() > 0) {
      *
      *       // To generate a star from a root wave rootItr,
      *       // loop over symmetry operations of space group.
      *       For each symmetry operation group[j] {
      *         Compute vec = (rootItr->indicesBz)*group[j]
      *         Set phase = rootItr->indicesBz .dot. (group[j].t)
      *         Check for cancellation of the star, set "cancel" flag
      *         Add wave to std::set<TWave> star if not added before
      *         // Here, use of a std::set simplifies test of uniqueness
      *       }
      *
      *       Copy all waves from star to std::vector<TWave> tempStar
      *       Sort tempStar by indicesBz, in descending order
      *       // Here, use of a std::vector for tempStar allows sorting
      *
      *       // Add waves in star to tempList and remove from list
      *       For each wave in tempStar {
      *         Append the wave to GArray<TWave> tempList
      *         Erase the wave from std::set<TWave> list
      *         // Here, use of a std::set for list simplifies erasure
      *       }
      *
      *       Initialize a Star object named newStar
      *       Assign values to members beginId, endId, size, cancel
      *
      *       // Assign values of newStar.invertFlag, rootItr, nextInvert
      *       if (nextInvert == -1) {
      *          // This is the second star in pair
      *          newStar.invertFlag = -1;
      *          nextInvert = 1;
      *          Set rootItr to the first wave in remaining list
      *       } else {
      *          Search for inverse of rootItr in this star
      *          if inverse is in this star {
      *             newStar.invertFlag = 0
      *             nextInvert = 1;
      *             Set rootItr to the first wave in remaining list
      *          } else
      *          Search for inverse of rootItr in this remaining list
      *          if the inverse is in the remaining list {
      *             newStar.invertFlag = 1
      *             nextInvert = -1;
      *             set rootItr to inverse of current root
      *          }
      *       }
      *
      *       Append newStar object to GArray<Star> stars_
      *
      *     } // end loop over stars in a single list
      *
      *     // At this point, tempList contains the contents of the
      *     // waves_ array occupying the range [beginId, endId-1],
      *     // grouped by stars, with waves within each star sorted
      *     // by indexBz.
      *
      *     // Overwrite the block of array waves_ with indices in the
      *     // range [beginId, endId-1] with the contents of tempList.
      *     For each wave in tempList {
      *       Copy a TWave in tempList to a Basis:Wave in  waves_
      *       Assign a complex coefficient of unit norm to the Wave
      *     }
      *
      *     // At this point, coefficients of waves have unit magnitude
      *     // and correct relative phases within each star, but not the
      *     // final absolute phases or magnitude.
      *
      *   } // finish processing of one list (waves of equal norm)
      *
      * } // End initial processing of all waves and stars
      *
      * // Set phases of wave coefficients
      * For each star in array stars_ {
      *   if star is closed under inversion (star.invertFlag == 0) {
      *     Find the inverse of every wave in this star
      *     if star is cancelled {
      *       Set coefficients of all waves to zero
      *     } else {
      *       Set the root to the first wave in the star
      *       For each wave in star:
      *          Divide coeff by the root coefficient
      *       }
      *       if (coeffs of root & inverse are not complex conjugates){
      *          Divide all coeffs by a common phasor chosen to obtain
      *             complex conjugate coefficients for root and inverse
      *       }
      *     }
      *   } else
      *   if (star.invertFlag == 1) {
      *     Find the inverse of every wave in this star and next star
      *     Set root of this star to the 1st wave in this star
      *     Set partner to the inverse of the root of this star
      *     If this star is cancelled {
      *       Set coefficients in this star and next to zero
      *     } else {
      *       For each wave in this star:
      *          Divide coeff by the root coefficient
      *       }
      *       For each wave in the next star:
      *          Divide coeff by the partner coefficient
      *       }
      *     }
      *   }
      *   // Note: If star.invertFlag = -1, do nothing because properties
      *   // of this star were all set when processing its partner.
      * }
      *
      * // For all waves, normalize coefficients and set starId
      * For each star in array stars_ {
      *   For each wave in this star {
      *     Set Wave.starId
      *     Divide coefficient by sqrt(double(star.size))
      *   }
      * }
      *
      * // For all waves, set implicit member and add to look up table
      * For each wave in array waves_ {
      *   Set Wave::implicit attribute
      *   Set waveIds_[rank] = i
      * }
      */

      // Loop over all waves (initial processing of waves)
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
         if (newList) {

            // Copy waves of equal norm into std::set "list"
            list.clear();
            tempList.clear();
            for (j = listBegin; j < listEnd; ++j) {
               wave.indicesDft = waves_[j].indicesDft;
               wave.indicesBz = waves_[j].indicesBz;
               wave.sqNorm = waves_[j].sqNorm;
               if (j > listBegin) {
                  UTIL_CHECK( std::abs(wave.sqNorm-waves_[j].sqNorm)
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
               rootVecDft = rootItr->indicesDft;
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
                  UTIL_CHECK(std::abs(Gsq - unitCell().ksq(vec)) < epsilon);

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

                  if (wave.indicesDft == rootVecDft) {
                     if (std::abs(wave.phase) > 1.0E-6) {
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
                     if (std::abs(phase_diff) > 1.0E-6) {
                        cancel = true;
                     }

                  }

               }

               // Copy all waves from std::set star to std::vector tempStar
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
               // basis functions (nBasis_) & waves in basis (nBasisWave_)
               if (!cancel) {
                  ++nBasis_;
                  nBasisWave_ += star.size();
               }

               // Initialize a Star object
               // newStar.eigen = Gsq;
               newStar.beginId = starBegin;
               newStar.endId = newStar.beginId + star.size();
               newStar.size = star.size();
               newStar.cancel = cancel;
               // Note: newStar.starInvert is not yet known

               // Determine invertFlag, rootItr and nextInvert
               if (nextInvert == -1) {

                  // If this star is 2nd of a pair related by inversion,
                  // set root for next star to 1st wave of remaining list.

                  newStar.invertFlag = -1;
                  rootItr = list.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the 2nd of a pair of partners,
                  // then determine if it is closed under inversion.

                  // Compute inverse nVec of root vector in FBZ
                  nVec.negate(rootVecBz);

                  // Shift inverse nVec to the DFT mesh
                  (*meshPtr_).shift(nVec);

                  // Search for inverse of root vector within this star
                  bool inverseFound = false;
                  setItr = star.begin();
                  for ( ; setItr != star.end(); ++setItr) {
                     if (nVec == setItr->indicesDft) {
                        inverseFound = true;
                        break;
                     }
                  }

                  if (inverseFound) {

                     // If this star is closed under inversion, the root
                     // of next star is the 1st vector of remaining list.

                     newStar.invertFlag = 0;
                     rootItr = list.begin();
                     nextInvert = 1;

                  } else {

                     // This star is open under inversion, and is the
                     // first star of a pair related by inversion

                     newStar.invertFlag = 1;
                     nextInvert = -1;

                     // Find inverse of the root of this star in the
                     // remaining list, and use this inverse as the
                     // root of the next star.

                     setItr = list.begin();
                     for ( ; setItr != list.end(); ++setItr) {
                        if (nVec == setItr->indicesDft) {
                           inverseFound = true;
                           rootItr = setItr;
                           break;
                        }
                     }
                     // If inverseFound, then rootVecDft = nVec

                     // Failure to find the inverse here is an error:
                     // It must be either in this star or remaining list

                     if (!inverseFound) {
                        std::cout << "Inverse not found for: " << "\n";
                        std::cout << " vec (ft):"
                                  << rootVecDft <<"\n";
                        std::cout << " vec (bz):"
                                  << rootVecBz <<"\n";
                        std::cout << "-vec (dft):" << nVec << "\n";
                        UTIL_CHECK(inverseFound);
                     }

                  }

               }

               stars_.append(newStar);
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
               if (std::abs(imag(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(real(coeff), 0.0);
               }
               if (std::abs(real(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(0.0, imag(coeff));
               }
               waves_[k].coeff = coeff;
            }

            // Processing of list is now complete.
            // Here, waves_[k].coeff has unit absolute magnitude, and
            // correct relative phases for waves within a star, but
            // the coeff may not be unity for the first or last wave in
            // the star.

            listBegin = listEnd;
         }
         // Finished processing a list of waves of equal norm

      } // End loop over all waves
      nStar_ = stars_.size();
      // Complete initial processing of all lists and stars

      /*
      * Conventions for phases of wave coefficients (imposed below):
      *   - Coefficients of the root of each star and its inverse must
      *     be complex conjugates.
      *   - In a closed star (starInvert = 0), the coefficient of the
      *     root must have a non-negative real part. If the root
      *     coefficient is pure imaginary, the imaginary part must
      *     be negative.
      *   - In a pair of open stars that are related by inversion
      *     symmetry, the coefficients of the first wave of the first 
      *     star and its inverse must have real coefficients.
      */

      // Final processing of phases of of waves in stars:
      std::complex<double> rootCoeff; // Coefficient of root wave
      std::complex<double> partCoeff; // Coefficient of partner of root
      std::complex<double> d;
      int rootId, partId;
      for (i = 0; i < nStar_; ++i) {

         // Treat open and closed stars differently

         if (stars_[i].invertFlag == 0) {

            // First, assign inverseId for each wave in star
            for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
               if (waves_[j].inverseId < 0) { // if inverseId is unassigned

                  // Compute nVec = inverse of root, shifted to DFT mesh
                  nVec.negate(waves_[j].indicesBz);
                  (*meshPtr_).shift(nVec);

                  // Find inverse

                  // Check in the position that the inverse is
                  // expected to be located for a typical star
                  k = stars_[i].endId - 1 - (j - stars_[i].beginId);
                  if (nVec == waves_[k].indicesDft) {
                     waves_[j].inverseId = k;
                     waves_[k].inverseId = j;
                  } else {
                     // Inverse not in expected position, search full star
                     // (this usually occurs for stars on the edge of the
                     // Brillouin zone)
                     for (k = j; k < stars_[i].endId; ++k) {
                        if (nVec == waves_[k].indicesDft) {
                           waves_[j].inverseId = k;
                           waves_[k].inverseId = j;
                           break;
                        }
                     }
                  }

                  // For invertFlag == 0, failure to find nVec in this
                  // star is a fatal error
                  if (waves_[j].inverseId < 0) {
                     std::cout << "\n";
                     std::cout << "Inverse not found in closed star"
                              << std::endl;
                     std::cout << "G = " << waves_[j].indicesBz
                               << ", coeff = " << waves_[j].coeff
                               << std::endl;
                     std::cout << "All waves in star " << i
                               << std::endl;
                     for (k=stars_[i].beginId; k < stars_[i].endId; ++k) {
                        std::cout << waves_[k].indicesBz << "  "
                                 << waves_[k].coeff << std::endl;
                     }
                     UTIL_CHECK(waves_[j].inverseId >= 0);
                  }

               }
            }

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

               // Set partId to index of the inverse of the root
               partId = waves_[rootId].inverseId;

               // Divide all coefficients by the root coefficient
               rootCoeff = waves_[rootId].coeff;
               for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
                  waves_[j].coeff /= rootCoeff;
               }
               rootCoeff = waves_[rootId].coeff;
               UTIL_CHECK(std::abs(real(rootCoeff) - 1.0) < 1.0E-9);
               UTIL_CHECK(std::abs(imag(rootCoeff)) < 1.0E-9);

               // Require coefficients of root and inverse are conjugates
               if (partId != rootId) {

                  partCoeff = waves_[partId].coeff;
                  UTIL_CHECK(std::abs(std::abs(partCoeff) - 1.0) < 1.0E-9);
                  if (std::abs(partCoeff - rootCoeff) > 1.0E-6) {
                     d = sqrt(partCoeff);
                     if (real(d) < -1.0E-4) {
                        d = -d;
                     } else
                     if (std::abs(real(d)) <= 1.0E-4) {
                        if (imag(d) < 0.0) {
                           d = -d;
                        }
                     }
                     for (j=stars_[i].beginId; j < stars_[i].endId; ++j){
                        waves_[j].coeff /= d;
                     }
                  }

               }

            } // end if (cancel) ... else ...

         } // end if (stars_[i].invertFlag == 0)
         else
         if (stars_[i].invertFlag == 1) {

            // Process a pair of open stars related by inversion.

            // Preconditions:
            UTIL_CHECK(stars_[i].size == stars_[i+1].size);
            UTIL_CHECK(stars_[i].cancel == stars_[i+1].cancel);
            // UTIL_CHECK(stars_[i+1].invertFlag == -1);

            // First, assign inverseId for each wave in pair of stars
            for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
               if (waves_[j].inverseId < 0) { // if inverseId is unassigned

                  // Compute nVec = inverse of root, shifted to DFT mesh
                  nVec.negate(waves_[j].indicesBz);
                  (*meshPtr_).shift(nVec);

                  // Find inverse

                  // Check in the position that the inverse is expected
                  // to be located for a typical pair of stars
                  k = stars_[i+1].endId - 1 - (j - stars_[i].beginId);
                  if (nVec == waves_[k].indicesDft) {
                     waves_[j].inverseId = k;
                     waves_[k].inverseId = j;
                  } else {
                     // Inverse not in expected position, search full star
                     // (this usually occurs for stars on the edge of the
                     // Brillouin zone)
                     k = stars_[i+1].beginId;
                     for ( ; k < stars_[i+1].endId; ++k) {
                        if (nVec == waves_[k].indicesDft) {
                           waves_[j].inverseId = k;
                           waves_[k].inverseId = j;
                           break;
                        }
                     }
                  }

                  // For invertFlag = 1, failure to find nVec in the
                  // next star is a fatal error
                  UTIL_CHECK(waves_[j].inverseId >= 0);

               }
            }

            // Check that inverseId was assigned for all waves in the
            // next star
            for (j = stars_[i+1].beginId; j < stars_[i+1].endId; ++j) {
               UTIL_CHECK(waves_[j].inverseId >= 0);
            }

            // Identify root of this star (star i)
            // Set the root to be the first wave in the star
            rootId = stars_[i].beginId;
            stars_[i].waveBz = waves_[rootId].indicesBz;

            // Identify root of the next star (star i+1)
            // Set the root to be the inverse of the root of star i
            partId = waves_[rootId].inverseId;
            stars_[i+1].waveBz = waves_[partId].indicesBz;

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

      // Set tiny real and imaginary parts to zero (due to round-off)
      for (i = 0; i < nWave_; ++i) {
         if (std::abs(real(waves_[i].coeff)) < 1.0E-8) {
            waves_[i].coeff
                     = std::complex<double>(0.0, imag(waves_[i].coeff));
         }
         if (std::abs(imag(waves_[i].coeff)) < 1.0E-8) {
            waves_[i].coeff
                    = std::complex<double>(real(waves_[i].coeff), 0.0);
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

      // Set Star::starId and Star::basisId, add to look-up table starIds_.
      starIds_.allocate(nBasis_);
      j = 0;
      for (i = 0; i < nStar_; ++i) {
         stars_[i].starId = i;
         if (stars_[i].cancel) {
            stars_[i].basisId = -1;
         } else {
            stars_[i].basisId = j;
            UTIL_CHECK(j < nBasis_);
            starIds_[j] = i;
            ++j;
         }
      }
      UTIL_CHECK(j == nBasis_);

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
            out << Int(i, 8);
            for (j = 0; j < D; ++j) {
               out << Int(waves_[i].indicesBz[j], 5);
            }
            out << Int(waves_[i].starId, 6);
            out << "  " << Dbl(waves_[i].coeff.real(), 15);
            out << "  " << Dbl(waves_[i].coeff.imag(), 15);
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
      int i, j;
      for (i = 0; i < nStar_; ++i) {
         if (outputAll || (!stars_[i].cancel)) {
            out << Int(stars_[i].basisId, 6);  // basisId
            out << Int(i, 6);                  // starId
            out << Int(stars_[i].size, 5)
                << Int(stars_[i].beginId, 8)
                << Int(stars_[i].endId, 8)
                << Int(stars_[i].invertFlag, 4);
            if (outputAll) {
               out << Int(stars_[i].cancel, 4);
            }
            for (j = 0; j < D; ++j) {
               out << Int(stars_[i].waveBz[j], 6);
            }
            out << std::endl;
         }
      }
   }

   template <int D>
   bool Basis<D>::isValid() const
   {
      double Gsq;
      IntVec<D> v;
      int is, ib, iw, iwp, j;

      // Check total number of waves == # of grid points
      if (nWave_ != mesh().size()) {
         std::cout << "nWave != size of mesh" << std::endl;
         return false;
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
         if (std::abs(Gsq - waves_[iw].sqNorm) > 1.0E-8) {
            std::cout << "\n";
            std::cout << "Incorrect sqNorm:" << "\n"
                      << "wave.indicesBz = " << "\n"
                      << "wave.sqNorm    = " << waves_[iw].sqNorm << "\n"
                      << "|v|^{2}        = " << Gsq << "\n";
            return false;
         }

         // Check that wave indicesBz is an image of indicesDft
         mesh().shift(v);
         if (v != waves_[iw].indicesDft) {
            std::cout << "\n";
            std::cout << "shift(indicesBz) != indicesDft" << std::endl;
            return false;
         }

         // Compare Wave::starId to Star::beginId and Star::endId
         is = waves_[iw].starId;
         if (iw < stars_[is].beginId) {
            std::cout << "\n";
            std::cout << "Wave::starId < Star::beginId" << std::endl;
            return false;
         }
         if (iw >= stars_[is].endId) {
            std::cout << "\n";
            std::cout << "Wave::starId >= Star::endId" << std::endl;
            return false;
         }

         // Check that inverseId has been assigned
         if (waves_[iw].inverseId < 0) {
            std::cout << "\n";
            std::cout << "Wave::inverseId not assigned\n";
            std::cout << "G = " << waves_[iw].indicesBz << std::endl;
            return false;
         }

         // Check that inverseId points to the correct wave
         v.negate(waves_[iw].indicesBz);
         mesh().shift(v);
         iwp = waves_[iw].inverseId;
         if (waves_[iwp].indicesDft != v) {
            std::cout << "\n";
            std::cout << "Wave::inverseId is not inverse" << std::endl;
            std::cout << "G = " << waves_[iw].indicesBz << std::endl;
            std::cout << "-G (from inverseId) = "
                      << waves_[iwp].indicesBz << std::endl;
            return false;
         }

         // Check that inverseId of the inverse wave is correct
         if (waves_[iwp].inverseId != iw) {
            std::cout << "\n";
            std::cout << "Wave::inverseId values do not agree\n";
            std::cout << "+G = " << waves_[iw].indicesBz << std::endl;
            std::cout << "-G = " << waves_[iwp].indicesBz << std::endl;
            return false;
         }

         // Check that either this wave or its inverse is explicit
         if (waves_[iw].implicit == true && waves_[iwp].implicit == true)
         {
            std::cout << "\n";
            std::cout << "Wave and its inverse are both implicit";
            std::cout << "+G = " << waves_[iw].indicesBz << std::endl;
            std::cout << "-G = " << waves_[iwp].indicesBz << std::endl;
            return false;
         }
      }

      // Loop over all stars (elements of stars_ array)
      int nWave = 0;
      for (is = 0; is < nStar_; ++is) {

         // Check star size
         nWave += stars_[is].size;
         if (stars_[is].size != stars_[is].endId - stars_[is].beginId) {
            std::cout << "\n";
            std::cout << "Inconsistent Star::size:" << std::endl;
            std::cout << "Star id    "  << is << std::endl;
            std::cout << "star size  "  << stars_[is].size << std::endl;
            std::cout << "Star begin "  << stars_[is].beginId << std::endl;
            std::cout << "Star end   "  << stars_[is].endId << std::endl;
            return false;
         }
         if (is > 0) {
            if (stars_[is].beginId != stars_[is-1].endId) {
               std::cout << "\n";
               std::cout << "Star ranges not consecutive:" << std::endl;
               std::cout << "Star id    "     << is << std::endl;
               std::cout << "stars_[" << is << "]"   << ".beginId = "
                         << stars_[is].beginId << std::endl;
               std::cout << "stars_[" << is-1 << "]" << ".endId   = "
                         << stars_[is-1].endId << std::endl;
               return false;
            }
         }

         // Check waveBz indices of star
         if (stars_[is].invertFlag == -1) {
            v.negate(stars_[is-1].waveBz);
            v = shiftToMinimum(v, mesh().dimensions(), *unitCellPtr_);
            if (stars_[is].waveBz != v) {
               std::cout << "\n";
               std::cout << "waveBz of star is not inverse of waveBz "
                         << "of previous star" << std::endl;
               std::cout << "star id " << is << std::endl;
               std::cout << "waveBz  " << stars_[is].waveBz << std::endl;
               std::cout << "waveBz (previous star) "
                         << stars_[is-1].waveBz << std::endl;
               return false;
            }
         } else {
            v = waves_[stars_[is].beginId].indicesBz;
            if (stars_[is].waveBz != v) {
               std::cout << "\n";
               std::cout << "waveBz of star != first wave of star"
                         << std::endl;
               std::cout << "star id    " << is << std::endl;
               std::cout << "waveBz     " << stars_[is].waveBz
                         << std::endl;
               std::cout << "first wave " << v << std::endl;
               return false;
            }
         }

         // Check star ids of waves in star
         for (iw = stars_[is].beginId; iw < stars_[is].endId; ++iw) {
            if (waves_[iw].starId != is) {
               std::cout << "\n";
               std::cout << "Inconsistent Wave::starId :" << std::endl;
               std::cout << "star id      "  << is << std::endl;
               std::cout << "star beginId "  << stars_[is].beginId << "\n";
               std::cout << "star endId   "  << stars_[is].endId << "\n";
               std::cout << "wave id      "  << iw << "\n";
               std::cout << "wave starId  "  << waves_[iw].starId << "\n";
               return false;
            }
         }

         // Check Star::starId is equal to array index
         if (stars_[is].starId != is) {
            std::cout << "\n";
            std::cout << "stars_[is].starId != is for "
                      << "is = " << is << "\n";
            return false;
         }

         // Check Star::basisId and starIds_ look up table
         ib = stars_[is].basisId;
         if (stars_[is].cancel) {
            if (ib != -1) {
               std::cout << "\n";
               std::cout << "basisId != -1 for cancelled star\n";
               std::cout << "star id = " << is << "\n";
               return false;
            }
         } else {
            if (starIds_[ib] != is) {
               std::cout << "\n";
               std::cout << "starIds_[stars_[is].basisId] != is for: \n";
               std::cout << "is                      = " << is << "\n";
               std::cout << "ib = stars_[is].basisId = " << ib << "\n";
               std::cout << "starIds_[ib]            = " << starIds_[ib]
                         << "\n";
               return false;
            }
         }

         // Check ordering of waves in star
         for (iw = stars_[is].beginId + 1; iw < stars_[is].endId; ++iw) {
            if (waves_[iw].indicesBz > waves_[iw-1].indicesBz) {
               std::cout << "\n";
               std::cout << "Failure of ordering by indicesB within star"
                         << std::endl;
               return false;
            }
            if (waves_[iw].indicesBz == waves_[iw-1].indicesBz) {
               std::cout << "\n";
               std::cout << "Equal values of indicesBz within star"
                         << std::endl;
               return false;
            }
         }

         // Check that all coefficients are zero if star is cancelled
         if (stars_[is].cancel) {
            for (iw = stars_[is].beginId + 1; iw < stars_[is].endId; ++iw)
            {
               if (std::abs(waves_[iw].coeff) > 1.0E-8) {
                  std::cout << "\n";
                  std::cout << "Nonzero coefficient in a cancelled star"
                              << "\n";
                  std::cout << "G = " << waves_[iw].indicesBz
                              << "  coeff = " << waves_[iw].coeff
                              << "\n";
                  return false;
               }
            }
         }

      } // End do loop over all stars

      // Check that all waves in mesh are accounted for in stars
      if (stars_[nStar_-1].endId != mesh().size()) {
         std::cout << "\n";
         std::cout << "Star endId of last star != mesh size" << std::endl;
         return false;
      }
      if (nWave != mesh().size()) {
         std::cout << "\n";
         std::cout << "Sum of star sizes != mesh size" << std::endl;
         return false;
      }

      // Loop over closed stars and related pairs of stars.
      // Test closure under inversion and conjugacy of coefficients.
      std::complex<double> cdel;
      bool cancel;
      is = 0;
      while (is < nStar_) {
         cancel = stars_[is].cancel;

         if (stars_[is].invertFlag == 0) {

            // Test that star is closed under inversion and real
            int begin = stars_[is].beginId;
            int end = stars_[is].endId;
            for (iw = begin; iw < end; ++iw) {
               iwp = waves_[iw].inverseId;
               if (waves_[iwp].starId != is) {
                  std::cout << "\n";
                  std::cout << "Inverse not found in closed star"
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
               if (!cancel) {
                  cdel = std::conj(waves_[iwp].coeff);
                  cdel -= waves_[iw].coeff;
                  if (std::abs(cdel) > 1.0E-8) {
                     std::cout << "\n";
                     std::cout << "Function for closed star is not real:"
                              << "\n";
                     std::cout << "+G = " << waves_[iw].indicesBz
                              << "  coeff = " << waves_[iw].coeff
                              << "\n";
                     std::cout << "-G = " << waves_[iwp].indicesBz
                              << "  coeff = " << waves_[iwp].coeff
                              << "\n";
                     std::cout << "Coefficients are not conjugates."
                               << "\n";
                     std::cout << "All waves in star " << is << "\n";
                     for (j=begin; j < end; ++j) {
                        std::cout << waves_[j].indicesBz << "  "
                                 << waves_[j].coeff << "\n";
                     }
                     return false;
                  }
               }
            }

            // Finished processing a closed star, increment counter is
            ++is;

         } else {

            // Test pairs of open stars

            if (stars_[is].invertFlag != 1) {
               std::cout << "\n";
               std::cout << "Expected invertFlag == 1" << std::endl;
               return false;
            }
            if (stars_[is+1].invertFlag != -1) {
               std::cout << "\n";
               std::cout << "Expected invertFlag == -1" << std::endl;
               return false;
            }
            if (stars_[is+1].size != stars_[is].size) {
               std::cout << "\n";
               std::cout << "Partner stars of different size" << std::endl;
               return false;
            }
            if (stars_[is+1].cancel != stars_[is].cancel) {
               std::cout << "\n";
               std::cout << "Partners stars with different cancel flags"
                         << std::endl;
               return false;
            }

            // Begin and end wave ids for the first and second stars
            int begin1 = stars_[is].beginId;
            int end1 = stars_[is].endId;
            int begin2 = stars_[is+1].beginId;
            int end2 = stars_[is+1].endId;

            // Check that inverse is in next star and check for
            // conjugate coefficients

            // Loop over waves in first star
            for (iw = begin1; iw < end1; ++iw) {
               iwp = waves_[iw].inverseId;
               if (waves_[iwp].starId != is + 1) {
                  std::cout << "\n";
                  std::cout << "Inverse not found for G in open star"
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
               }
               if (!cancel) {
                  cdel = std::conj(waves_[iwp].coeff);
                  cdel -= waves_[iw].coeff;
                  if (std::abs(cdel) > 1.0E-8) {
                     std::cout << "\n";
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
            }

            // Finished processing a pair, increment star counter by 2
            is += 2;

         } // end if (stars_[is].invertFlag == 0) ... else ...

      } // end while (is < nStar_) loop over stars

      // Loop over basis functions
      for (ib = 0; ib < nBasis_; ++ib) {
         is = starIds_[ib];
         if (stars_[is].cancel) {
            std::cout << "\n";
            std::cout << "Star referred to by starIds_ is cancelled\n";
            return false;
         }
         if (stars_[is].basisId != ib) {
            std::cout << "\n";
            std::cout << "Error: stars_[starIds_[ib]].basisId != ib\n";
            std::cout << "Basis function index ib = " << ib << "\n";
            std::cout << "is = starIds_[ib]       = " << is << "\n";
            std::cout << "stars_[is].basisId      = "
                      << stars_[is].basisId << "\n";
            return false;
         }
      }

      // The end of this function is reached iff all tests passed.
      return true;
   }

   /*
   * Get the associated signal (notifies observers of basis initialization).
   */
   template <int D>
   Signal<void>& Basis<D>::signal()
   {
      UTIL_CHECK(signalPtr_); 
      return *signalPtr_;
   }

}
}
#endif
