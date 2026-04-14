#ifndef PRDC_BASIS_TPP
#define PRDC_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"
#include "BWave.h"
#include "groupFile.h"
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/signal/Signal.h>
#include <util/format/Dbl.h>
#include <util/math/Constants.h>

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
   *  - Member array waves_ contains list of waves ordered by sqNorm.
   *  - Each wave has indicesStd, indicesMin and sqNorm set.
   *  - Member array stars_  is still empty.
   */
   template <int D>
   void Basis<D>::makeWaves()
   {
      IntVec<D> meshDimensions = mesh().dimensions();
      std::vector< BWave<D> > twaves;
      twaves.reserve(nWave_);

      // Loop over k-grid mesh to generate all waves, add to twaves
      BWave<D> w;
      IntVec<D> v;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         w.indicesStd = itr.position();
         v = shiftToMinimum(w.indicesStd, meshDimensions, *unitCellPtr_);
         w.indicesMin = v;
         w.sqNorm = unitCell().ksq(v);
         twaves.push_back(w);
      }

      // Sort temporary array twaves
      BWaveNormComp<D> comp;
      std::sort(twaves.begin(), twaves.end(), comp);

      // Copy temporary array twaves into member variable waves_
      for (int i = 0; i < nWave_; ++i) {
         waves_[i].sqNorm = twaves[i].sqNorm;
         waves_[i].indicesStd = twaves[i].indicesStd;
         waves_[i].indicesMin = twaves[i].indicesMin;
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
      *  - A "bunch" is a set of wavevectors of equal magnitude.
      *  - A "star" is a set of wavevectors that are related by symmetry.
      *
      * Each bunch contains one or more complete stars. Bunches are
      * identified as an intermediate step in identification of stars.
      *
      * During initial processing, wavevectors are temporarily stored in 
      * BWave<D> objects.  The following local containers of BWave<D> 
      * objects are used:
      *
      *   bunch - a set of waves of equal norm 
      *   star - a set of symmetry-related waves
      *   starList - a sorted star, sorted by descending indicesMin
      *   bunchList - a sorted bunch, with contiguous sorted stars
      *
      * The bunch and star containers are std::set< BWave<D> > objects. 
      * The starList container is a std::vector< BWave<D> > objects. 
      * The use of std::set simplifies addition of waves, by allowing
      * easy identification of duplicates, as well as later removal.
      * The use of std::vector for the starList allows the use of the 
      * C++ std::sort algorithm to sort stars by standard wavevector 
      * indices.
      */

      // Local BWave<D> containers and associated iterators
      std::set< BWave<D>, BWaveStdComp<D> > bunch;
      std::set< BWave<D>, BWaveStdComp<D> > star;
      std::vector< BWave<D> > starList;
      GArray< BWave<D> > bunchList;
      typename std::set< BWave<D>, BWaveStdComp<D> >::iterator rootItr;
      typename std::set< BWave<D>, BWaveStdComp<D> >::iterator setItr;

      // Local variables
      BWave<D> wave;
      Basis<D>::Star newStar;
      std::complex<double> coeff;
      double Gsq;
      double Gsq_max;
      double phase_diff;
      const double twoPi = 2.0*Constants::Pi;
      const double epsilon = 1.0E-8;
      IntVec<D> meshDimensions = mesh().dimensions();
      IntVec<D> rootVecMin;  // Min indices for root of this star
      IntVec<D> rootVecStd;  // Std indices for root of this star
      IntVec<D> vec;         // Indices of temporary wavevector
      IntVec<D> nVec;        // Indices of inverse of a wavevector
      int bunchBegin = 0;    // id of first wave in this bunch 
      int bunchEnd = 0;      // (id of last wave in this bunch) + 1
      int bunchSize;         // bunchEnd - bunchBegin
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
      *   Search for end of a "bunch" (i.e., contiguous block of waves
      *   of equal magnitude) by identifying changes in magnitude. The
      *   resulting bunch has indices bunchBegin, ... , bunchEnd-1 .
      *   Set newList true.
      *
      *   // Each bunch may contain one or more stars.
      *   // Process the newly identified bunch to identify stars.
      *
      *   If (newList) {
      *
      *     Copy all waves in the range into std::set bunch
      *
      *     Set rootItr to the first wave in the bunch
      *
      *     // Loop over stars within the bunch
      *     while (bunch.size() > 0) {
      *
      *       // To generate a star from a root wave rootItr,
      *       // loop over symmetry operations of space group.
      *
      *       For each symmetry operation group[j] {
      *         Compute vec = (rootItr->indicesMin)*group[j]
      *         Set phase = rootItr->indicesMin .dot. (group[j].t)
      *         Check for cancellation of the star, set "cancel" flag
      *         Add wave to std::set<BWave> star if not added before
      *         // Here, use of a std::set simplifies test of uniqueness
      *       }
      *
      *       Copy all waves from star to std::vector<BWave> starList
      *       Sort starList by indicesMin, in descending order
      *       // Here, use of a std::vector for starList allows sorting
      *
      *       // Add waves in star to bunchList and remove from bunch
      *       For each wave in starList {
      *         Append the wave to GArray<BWave<D> > bunchList
      *         Erase the wave from std::set<BWave<D> > bunch
      *         // Here, use of a std::set for bunch simplifies erasure
      *       }
      *
      *       Initialize a Star object named newStar
      *       Assign values to members beginId, endId, size, cancel
      *
      *       // Assign values of newStar.invertFlag, nextInvert, rootItr
      *       if (nextInvert == -1) {
      *          // This is the second star in pair
      *          newStar.invertFlag = -1;
      *          nextInvert = 1;
      *          Set rootItr to the first wave in remaining bunch
      *       } else {
      *          Search for inverse of rootItr in this star
      *          if inverse is in this star {
      *             // This is a closed star
      *             newStar.invertFlag = 0
      *             nextInvert = 1;
      *             Set rootItr to the first wave in remaining bunch
      *          } else
      *          Search for inverse of rootItr in this remaining bunch
      *          if the inverse is in the remaining bunch {
      *             // This is the first open star in a pair
      *             newStar.invertFlag = 1
      *             nextInvert = -1;
      *             set rootItr to inverse of current root
      *          }
      *       }
      *
      *       Append newStar object to GArray<Star> stars_
      *
      *     } // end loop over stars in a single bunch
      *
      *     // At this point, bunchList contains the contents of the
      *     // waves_ array occupying the range [beginId, endId-1],
      *     // grouped by stars, with waves within each star sorted
      *     // by minimal indices.
      *
      *     // Overwrite the block of array waves_ with indices in the
      *     // range [beginId, endId-1] with the contents of bunchList.
      *     For each wave in bunchList {
      *        Copy a BWave in bunchList to a Basis:Wave in  waves_
      *        Assign a complex coefficient of unit norm to the Wave
      *     }
      *
      *     // At this point, coefficients of waves have unit magnitude
      *     // and correct relative phases within each star, but not the
      *     // final absolute phases or magnitude.
      *
      *   } // finish processing of one bunch (waves of equal norm)
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

         // Determine if this wave begins a new bunch
         bool newList = false;
         if (i == nWave_) {
            bunchEnd = i;
            bunchSize = bunchEnd - bunchBegin;
            newList = true;
         } else {
            Gsq = waves_[i].sqNorm;
            if (Gsq > Gsq_max + epsilon) {
               Gsq_max = Gsq;
               bunchEnd = i;
               bunchSize = bunchEnd - bunchBegin;
               newList = true;
            }
         }

         // Process completed bunch of wavectors of equal norm
         if (newList) {

            // Copy waves of equal norm into std::set "bunch"
            bunch.clear();
            bunchList.clear();
            for (j = bunchBegin; j < bunchEnd; ++j) {
               wave.indicesStd = waves_[j].indicesStd;
               wave.indicesMin = waves_[j].indicesMin;
               wave.sqNorm = waves_[j].sqNorm;
               if (j > bunchBegin) {
                  UTIL_CHECK( std::abs(wave.sqNorm-waves_[j].sqNorm)
                                 < 2.0*epsilon );
               }
               bunch.insert(wave);
            }

            // On entry to each iteration of the loop over stars,
            // rootItr and nextInvert are known. The iterator rootItr
            // points to the wave in the remaining bunch that will be
            // used as the root of the next star. The flag nextInvert
            // is equal to -1 iff the previous star was the first of
            // a pair that are open under inversion, and is equal
            // to + 1 otherwise.

            // Initial values for first star in this bunch
            rootItr = bunch.begin();
            int nextInvert = 1;

            // Loop over stars with a bunch of waves of equal norm,
            // removing each star from set bunch as it is identified.
            // The root of the next star must have been chosen on
            // entry to each iteration of this loop.

            while (bunch.size() > 0) {

               rootVecMin = rootItr->indicesMin;
               rootVecStd = rootItr->indicesStd;
               Gsq = rootItr->sqNorm;
               cancel = false;
               star.clear();

               // Construct a star from root vector, by applying every
               // symmetry operation in the group to the root wavevector.
               for (j = 0; j < group.size(); ++j) {

                  // Apply symmetry (i.e., multiply by rotation matrix)
                  // vec = rotated wavevector.
                  vec = rootVecMin*group[j];

                  // Check that rotated vector has same norm as root.
                  UTIL_CHECK(std::abs(Gsq - unitCell().ksq(vec)) < epsilon);

                  // Initialize BWave object associated with rotated wave
                  wave.sqNorm = Gsq;
                  wave.indicesMin = shiftToMinimum(vec, meshDimensions,
                                                  *unitCellPtr_);
                  wave.indicesStd = vec;
                  mesh().shift(wave.indicesStd);

                  // Compute phase for coeff. of wave in basis function.
                  // Convention -pi < phase <= pi.
                  wave.phase = 0.0;
                  for (k = 0; k < D; ++k) {
                     wave.phase += rootVecMin[k]*(group[j].t(k));
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

                  if (wave.indicesStd == rootVecStd) {
                     if (std::abs(wave.phase) > 1.0E-6) {
                        cancel = true;
                     }
                  }

                  // Search for an equivalent wave already in the star.
                  // Note: Equivalent waves have equal standard indices.
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

               // Copy waves from std::set star to std::vector starList
               starList.clear();
               setItr = star.begin();
               for ( ; setItr != star.end(); ++setItr) {
                  starList.push_back(*setItr);
               }

               // Sort starList, in descending order by indicesMin.
               BWaveMinComp<D> waveMinComp;
               std::sort(starList.begin(), starList.end(), waveMinComp);

               // Append starList to bunchList, erase from set bunch
               int starListSize = starList.size();
               for (j = 0; j < starListSize; ++j) {
                  bunch.erase(starList[j]);
                  bunchList.append(starList[j]);
               }
               UTIL_CHECK(int(bunchList.size()+bunch.size())==bunchSize);

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
                  // set root for next star to 1st wave of remaining 
                  // bunch.

                  newStar.invertFlag = -1;
                  rootItr = bunch.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the 2nd of a pair of partners,
                  // then determine if it is closed under inversion.

                  // Compute negation nVec of root vector
                  nVec.negate(rootVecMin);

                  // Shift nVec to the standard mesh
                  (*meshPtr_).shift(nVec);

                  // Search for inverse of root vector within this star
                  bool inverseFound = false;
                  setItr = star.begin();
                  for ( ; setItr != star.end(); ++setItr) {
                     if (nVec == setItr->indicesStd) {
                        inverseFound = true;
                        break;
                     }
                  }

                  if (inverseFound) {

                     // If this star is closed under inversion, the root 
                     // of the next star is the 1st vector of remaining 
                     // bunch.

                     newStar.invertFlag = 0;
                     rootItr = bunch.begin();
                     nextInvert = 1;

                  } else {

                     // This star is open under inversion, and is the
                     // first star of a pair related by inversion

                     newStar.invertFlag = 1;
                     nextInvert = -1;

                     // Find inverse of the root of this star in the
                     // remaining bunch, and use this inverse as the
                     // root of the next star.

                     setItr = bunch.begin();
                     for ( ; setItr != bunch.end(); ++setItr) {
                        if (nVec == setItr->indicesStd) {
                           inverseFound = true;
                           rootItr = setItr;
                           break;
                        }
                     }
                     // If inverseFound, then rootVecStd = nVec

                     // Failure to find the inverse here is an error:
                     // It must be either in this star or remaining bunch

                     if (!inverseFound) {
                        std::cout << "Inverse not found for: " << "\n";
                        std::cout << " vec (std):"
                                  << rootVecStd <<"\n";
                        std::cout << " vec (min):"
                                  << rootVecMin <<"\n";
                        std::cout << "-vec (std):" << nVec << "\n";
                        UTIL_CHECK(inverseFound);
                     }

                  }

               }

               stars_.append(newStar);
               starBegin = newStar.endId;

            }
            // End loop over stars within a bunch.

            UTIL_CHECK(bunch.size() == 0);
            UTIL_CHECK(bunchList.size() == bunchEnd - bunchBegin);

            // Copy bunchList into corresponding section of waves_,
            // overwriting the section of waves_ used to create the bunch.
            // Compute a complex coefficient of unit norm for each wave.
            for (j = 0; j < bunchList.size(); ++j) {
               k = j + bunchBegin;
               waves_[k].indicesStd = bunchList[j].indicesStd;
               waves_[k].indicesMin = bunchList[j].indicesMin;
               waves_[k].sqNorm = bunchList[j].sqNorm;
               coeff = std::complex<double>(0.0, bunchList[j].phase);
               coeff = exp(coeff);
               if (std::abs(imag(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(real(coeff), 0.0);
               }
               if (std::abs(real(coeff)) < 1.0E-6) {
                  coeff = std::complex<double>(0.0, imag(coeff));
               }
               waves_[k].coeff = coeff;
            }

            // Processing of this bunch is now complete.
            // Here, waves_[k].coeff has unit absolute magnitude, and
            // correct relative phases for waves within a star, but 
            // the coeff may not be unity for the first or last wave
            // in the star.

            bunchBegin = bunchEnd;
         }
         // Finished processing a bunch of waves of equal norm

      } // End loop over all waves
      nStar_ = stars_.size();
      // Complete initial processing of all bunches and stars

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

                  // Compute nVec = inverse of root, shifted to std mesh
                  nVec.negate(waves_[j].indicesMin);
                  (*meshPtr_).shift(nVec);

                  // Find inverse

                  // Check in the position that the inverse is
                  // expected to be located for a typical star
                  k = stars_[i].endId - 1 - (j - stars_[i].beginId);
                  if (nVec == waves_[k].indicesStd) {
                     waves_[j].inverseId = k;
                     waves_[k].inverseId = j;
                  } else {
                     // Inverse not in expected position, search full star
                     // (this usually occurs for stars on the edge of the
                     // Brillouin zone)
                     for (k = j; k < stars_[i].endId; ++k) {
                        if (nVec == waves_[k].indicesStd) {
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
                     std::cout << "G = " << waves_[j].indicesMin
                               << ", coeff = " << waves_[j].coeff
                               << std::endl;
                     std::cout << "All waves in star " << i
                               << std::endl;
                     for (k=stars_[i].beginId; k < stars_[i].endId; ++k) {
                        std::cout << waves_[k].indicesMin << "  "
                                 << waves_[k].coeff << std::endl;
                     }
                     UTIL_CHECK(waves_[j].inverseId >= 0);
                  }

               }
            }

            // Identify root of this star (star i)
            // Set the root to be the first wave in the star

            rootId = stars_[i].beginId;
            stars_[i].waveMin = waves_[rootId].indicesMin;

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

                  // Compute nVec = inverse of root, shifted to std mesh
                  nVec.negate(waves_[j].indicesMin);
                  (*meshPtr_).shift(nVec);

                  // Find inverse

                  // Check in the position that the inverse is expected
                  // to be located for a typical pair of stars
                  k = stars_[i+1].endId - 1 - (j - stars_[i].beginId);
                  if (nVec == waves_[k].indicesStd) {
                     waves_[j].inverseId = k;
                     waves_[k].inverseId = j;
                  } else {
                     // Inverse not in expected position, search full star
                     // (this usually occurs for stars on the edge of the
                     // Brillouin zone)
                     k = stars_[i+1].beginId;
                     for ( ; k < stars_[i+1].endId; ++k) {
                        if (nVec == waves_[k].indicesStd) {
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
            stars_[i].waveMin = waves_[rootId].indicesMin;

            // Identify root of the next star (star i+1)
            // Set the root to be the inverse of the root of star i
            partId = waves_[rootId].inverseId;
            stars_[i+1].waveMin = waves_[partId].indicesMin;

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
         vec = waves_[i].indicesStd;

         // Validity check - check that vec is in standard k-grid mesh
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
               out << Int(waves_[i].indicesMin[j], 5);
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
               out << Int(stars_[i].waveMin[j], 6);
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

      // Loop over k-grid mesh to check consistency of waveIds_ and waves_
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         v = itr.position();
         iw = waveId(v);
         if (wave(iw).indicesStd != v) {
            std::cout << "Inconsistent waveId and Wave::indicesStd"
                      << std::endl;
            return false;
         }
      }

      // Loop over elements of waves_, check consistency of wave data.
      for (iw = 0; iw < nWave_; ++iw) {

         // Check sqNorm
         v = waves_[iw].indicesMin;
         Gsq = unitCell().ksq(v);
         if (std::abs(Gsq - waves_[iw].sqNorm) > 1.0E-8) {
            std::cout << "\n";
            std::cout << "Incorrect sqNorm:" << "\n"
                      << "wave.indicesMin = " << "\n"
                      << "wave.sqNorm    = " << waves_[iw].sqNorm << "\n"
                      << "|v|^{2}        = " << Gsq << "\n";
            return false;
         }

         // Check that wave indicesMin is an image of indicesStd
         mesh().shift(v);
         if (v != waves_[iw].indicesStd) {
            std::cout << "\n";
            std::cout << "shift(indicesMin) != indicesStd" << std::endl;
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
            std::cout << "G = " << waves_[iw].indicesMin << std::endl;
            return false;
         }

         // Check that inverseId points to the correct wave
         v.negate(waves_[iw].indicesMin);
         mesh().shift(v);
         iwp = waves_[iw].inverseId;
         if (waves_[iwp].indicesStd != v) {
            std::cout << "\n";
            std::cout << "Wave::inverseId is not inverse" << std::endl;
            std::cout << "G = " << waves_[iw].indicesMin << std::endl;
            std::cout << "-G (from inverseId) = "
                      << waves_[iwp].indicesMin << std::endl;
            return false;
         }

         // Check that inverseId of the inverse wave is correct
         if (waves_[iwp].inverseId != iw) {
            std::cout << "\n";
            std::cout << "Wave::inverseId values do not agree\n";
            std::cout << "+G = " << waves_[iw].indicesMin << std::endl;
            std::cout << "-G = " << waves_[iwp].indicesMin << std::endl;
            return false;
         }

         // Check that either this wave or its inverse is explicit
         if (waves_[iw].implicit == true && waves_[iwp].implicit == true)
         {
            std::cout << "\n";
            std::cout << "Wave and its inverse are both implicit";
            std::cout << "+G = " << waves_[iw].indicesMin << std::endl;
            std::cout << "-G = " << waves_[iwp].indicesMin << std::endl;
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

         // Check waveMin indices of star
         if (stars_[is].invertFlag == -1) {
            v.negate(stars_[is-1].waveMin);
            v = shiftToMinimum(v, mesh().dimensions(), *unitCellPtr_);
            if (stars_[is].waveMin != v) {
               std::cout << "\n";
               std::cout << "waveMin of star is not inverse of waveMin "
                         << "of previous star" << std::endl;
               std::cout << "star id " << is << std::endl;
               std::cout << "waveMin  " << stars_[is].waveMin << std::endl;
               std::cout << "waveMin (previous star) "
                         << stars_[is-1].waveMin << std::endl;
               return false;
            }
         } else {
            v = waves_[stars_[is].beginId].indicesMin;
            if (stars_[is].waveMin != v) {
               std::cout << "\n";
               std::cout << "waveMin of star != first wave of star"
                         << std::endl;
               std::cout << "star id    " << is << std::endl;
               std::cout << "waveMin    " << stars_[is].waveMin
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
            if (waves_[iw].indicesMin > waves_[iw-1].indicesMin) {
               std::cout << "\n";
               std::cout << "Failure of ordering by indicesB within star"
                         << std::endl;
               return false;
            }
            if (waves_[iw].indicesMin == waves_[iw-1].indicesMin) {
               std::cout << "\n";
               std::cout << "Equal values of indicesMin within star"
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
                  std::cout << "G = " << waves_[iw].indicesMin
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
                  std::cout << "G = " << waves_[iw].indicesMin
                            << "coeff = " << waves_[iw].coeff
                            << std::endl;
                  std::cout << "All waves in star " << is << "\n";
                  for (j=begin; j < end; ++j) {
                     std::cout << waves_[j].indicesMin << "  "
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
                     std::cout << "+G = " << waves_[iw].indicesMin
                              << "  coeff = " << waves_[iw].coeff
                              << "\n";
                     std::cout << "-G = " << waves_[iwp].indicesMin
                              << "  coeff = " << waves_[iwp].coeff
                              << "\n";
                     std::cout << "Coefficients are not conjugates."
                               << "\n";
                     std::cout << "All waves in star " << is << "\n";
                     for (j=begin; j < end; ++j) {
                        std::cout << waves_[j].indicesMin << "  "
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
                  std::cout << "+G = " << waves_[iw].indicesMin
                            << "coeff = " << waves_[iw].coeff
                            << std::endl;
                  std::cout << "Waves in star " << is
                            << "  (starInvert ==1):" << "\n";
                  for (j = begin1; j < end1; ++j) {
                     std::cout << waves_[j].indicesMin  << "  "
                               << waves_[j].coeff << "\n";
                  }
                  std::cout << "Waves in star " << is+1
                            << "  (starInvert == -1):" << "\n";
                  for (j=begin2; j < end2; ++j) {
                     std::cout << waves_[j].indicesMin  << "  "
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
                     std::cout << "+G = " << waves_[iw].indicesMin
                              << "  coeff = " << waves_[iw].coeff
                              << "\n";
                     std::cout << "-G = " << waves_[iwp].indicesMin
                              << "  coeff = " << waves_[iwp].coeff
                              << "\n";
                     std::cout << "Coefficients are not conjugates."
                               << "\n";
                     std::cout << "Waves in star " << is
                               << "  (starInvert ==1):" << "\n";
                     for (j = begin1; j < end1; ++j) {
                        std::cout << waves_[j].indicesMin  << "  "
                                  << waves_[j].coeff << "\n";
                     }
                     std::cout << "Waves in star " << is+1
                               << "  (starInvert == -1):" << "\n";
                     for (j=begin2; j < end2; ++j) {
                        std::cout << waves_[j].indicesMin  << "  "
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
