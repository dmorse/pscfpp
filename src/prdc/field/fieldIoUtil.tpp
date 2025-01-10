#ifndef PRDC_FIELD_IO_UTIL_TPP
#define PRDC_FIELD_IO_UTIL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2024, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "fieldIoUtil.h"

#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>

#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Prdc {

   /*
   * Check allocation of a single field of type FT, allocate if needed.
   */
   template <int D, class FT>
   void checkAllocateField(FT& field, 
                           IntVec<D> const& dimensions)
   {
      if (field.isAllocated()) {
         UTIL_CHECK(field.meshDimensions() == dimensions);
      } else {
         field.allocate(dimensions);
      }
   }

   /*
   * Check allocation of an array of fields of type FT, allocate if needed.
   */
   template <int D, class FT>
   void checkAllocateFields(DArray<FT>& fields,
                            IntVec<D> const& dimensions,
                            int nMonomer)
   {
      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();
         UTIL_CHECK(nMonomerFields > 0)
         UTIL_CHECK(nMonomerFields == nMonomer)
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].meshDimensions() == dimensions);
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(dimensions);
         }
      }
   }

   /*
   * Check allocation of an array of arrays of type AT, allocate if needed.
   */
   template <int D, class AT>
   void checkAllocateArrays(DArray<AT>& arrays,
                            int capacity,
                            int nMonomer)
   {
      if (arrays.isAllocated()) {
         int nMonomerArrays = arrays.capacity();
         UTIL_CHECK(nMonomerArrays > 0)
         UTIL_CHECK(nMonomerArrays == nMonomer)
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(arrays[i].capacity() == capacity);
         }
      } else {
         arrays.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            arrays[i].allocate(capacity);
         }
      }
   }

   // RGrid File Io Templates

   template <int D, class AT>
   void readRGridData(std::istream& in,
                      DArray<AT>& fields,
                      IntVec<D> const& dimensions,
                      int nMonomer)
   {
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(fields.capacity() == nMonomer);

      MeshIteratorFortran<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         for (int k = 0; k < nMonomer; ++k) {
            in >> fields[k][rank];
         }
      }
   }

   template <int D, class AT>
   void readRGridData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dimensions)
   {
      MeshIteratorFortran<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> field[rank];
      }
   }

   template <int D, class AT>
   void writeRGridData(std::ostream& out,
                       DArray<AT> const& fields,
                       IntVec<D> const& dimensions,
                       int nMonomer)
   {
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());

      MeshIteratorFortran<D> iter(dimensions);
      int rank, j;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         for (j = 0; j < nMonomer; ++j) {
            out << "  " << Dbl(fields[j][rank], 21, 13);
         }
         out << std::endl;
      }

   }

   template <int D, class AT>
   void writeRGridData(std::ostream& out,
                       AT const& field,
                       IntVec<D> const& dimensions)
   {
      MeshIteratorFortran<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         out << "  " << Dbl(field[rank], 21, 13);
         out << std::endl;
      }
   }

   // KGrid file IO templates

   template <int D, class AT>
   void readKGridData(std::istream& in,
                      DArray<AT> & fields,
                      IntVec<D> const& dftDimensions,
                      int nMonomer)
   {
      MeshIterator<D> iter(dftDimensions);
      int rank, i, j, idum;
      i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == rank);
         for (j = 0; j < nMonomer; ++j) {
            in >> fields[j][rank][0];
            in >> fields[j][rank][1];
         }
         ++i;
      }
   }

   template <int D, class AT>
   void readKGridData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dftDimensions)
   {
      MeshIterator<D> iter(dftDimensions);
      int rank, idum;
      int i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == rank);
         in >> field[rank][0];
         in >> field[rank][1];
         ++i;
      }
   }

   template <int D, class AT>
   void writeKGridData(std::ostream& out,
                       DArray<AT> const& fields,
                       IntVec<D> const& dftDimensions,
                       int nMonomer)
   {
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());

      MeshIterator<D> iter(dftDimensions);
      int rank;
      int i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         UTIL_CHECK(i == rank);
         out << Int(rank, 5);
         for (int j = 0; j < nMonomer; ++j) {
            out << "  "
                << Dbl(fields[j][rank][0], 21, 13)
                << Dbl(fields[j][rank][1], 21, 13);
         }
         out << std::endl;
         ++i;
      }

   }

   template <int D, class AT>
   void writeKGridData(std::ostream& out, 
                       AT const& field,
                       IntVec<D> const& dftDimensions)
   {
      MeshIterator<D> iter(dftDimensions);
      int rank, i;
      i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         UTIL_CHECK(i == rank);
         out << Int(rank, 5);
         out << "  "
             << Dbl(field[rank][0], 21, 13)
             << Dbl(field[rank][1], 21, 13);
         out << std::endl;
         ++i;
      }
   }

   /*
   * Read a set of fields in basis format.
   */
   template <int D>
   void readBasisData(std::istream& in,
                      DArray< DArray<double> >& fields,
                      UnitCell<D> const& unitCell,
                      Mesh<D> const& mesh,
                      Basis<D> const& basis,
                      int nStarIn)
   {
      UTIL_CHECK(basis.isInitialized());
      UTIL_CHECK(fields.isAllocated());

      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Initialize all field array elements to zero
      int fieldCapacity = fields[0].capacity();
      UTIL_CHECK(fieldCapacity > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].capacity() == fieldCapacity);
         for (int j = 0; j < fieldCapacity; ++j) {
            fields[i][j] = 0.0;
         }
      }

      // Reset nStarIn = min(nStarIn, fieldCapacity)
      if (fieldCapacity < nStarIn) {
         nStarIn = fieldCapacity;
      }

      // Allocate temp arrays used to read in components
      DArray<double> temp, temp2;
      temp.allocate(nMonomer);
      temp2.allocate(nMonomer);

      typename Basis<D>::Star const * starPtr;
      typename Basis<D>::Star const * starPtr2;
      IntVec<D> waveIn, waveIn2;
      int sizeIn, sizeIn2;
      int starId, starId2;
      int basisId, basisId2;
      int waveId, waveId2;

      std::complex<double> coeff, phasor;
      IntVec<D> waveBz, waveDft;
      int nReversedPair = 0;
      bool waveExists, sizeMatches;

      // Loop over stars in input file to read field components
      int i = 0;
      while (i < nStarIn) {

         // Read next line of data
         for (int j = 0; j < nMonomer; ++j) {
            in >> temp[j];        // field components
         }
         in >> waveIn;            // wave of star
         in >> sizeIn;            // # of waves in star
         ++i;

         sizeMatches = false;
         waveExists = false;

         // Check if waveIn is in first Brillouin zone (FBZ) for the mesh.
         waveBz = shiftToMinimum(waveIn, mesh.dimensions(), unitCell);
         waveExists = (waveIn == waveBz);

         if (waveExists) {

            // Find the star containing waveIn
            waveDft = waveIn;
            mesh.shift(waveDft);
            waveId = basis.waveId(waveDft);
            starId = basis.wave(waveId).starId;
            starPtr = &basis.star(starId);
            UTIL_CHECK(!(starPtr->cancel));
            basisId = starPtr->basisId;

            if (starPtr->size == sizeIn) {
               sizeMatches = true;
            } else {
               Log::file()
                  <<  "Warning: Inconsistent star size (line ignored)\n"
                  <<  "wave from file = " << waveIn << "\n"
                  <<  "size from file = " << sizeIn << "\n"
                  <<  "size of star   = " << starPtr->size
                  << "\n";
               sizeMatches = false;
            }

         }

         if (waveExists && sizeMatches) { // Attempt to process wave

            if (starPtr->invertFlag == 0) {

               if (starPtr->waveBz == waveIn) {

                  // Copy components of closed star to fields array
                  for (int j = 0; j < nMonomer; ++j) {
                      fields[j][basisId] = temp[j];
                  }

               } else {

                  Log::file()
                     <<  "Inconsistent wave of closed star on input\n"
                     <<  "wave from file = " << waveIn  << "\n"
                     <<  "starId of wave = " << starId  << "\n"
                     <<  "waveBz of star = " << starPtr->waveBz
                     << "\n";

               }

            } else {

               // Read the next line
               for (int j = 0; j < nMonomer; ++j) {
                  in >> temp2[j];          // components of field
               }
               in >> waveIn2;              // wave of star
               in >> sizeIn2;              // # of wavevectors in star
               ++i;

               // Check that waveIn2 is also in the 1st BZ
               waveBz =
                  shiftToMinimum(waveIn2, mesh.dimensions(), unitCell);
               UTIL_CHECK(waveIn2 == waveBz);

               // Identify the star containing waveIn2
               waveDft = waveIn2;
               mesh.shift(waveDft);
               waveId2 = basis.waveId(waveDft);
               starId2 = basis.wave(waveId2).starId;
               starPtr2 = &basis.star(starId2);
               UTIL_CHECK(!(starPtr2->cancel));
               basisId2 = starPtr2->basisId;
               UTIL_CHECK(starPtr2->size == sizeIn2);
               UTIL_CHECK(sizeIn == sizeIn2);

               if (starPtr->invertFlag == 1) {

                  // This is a pair of open stars written in the same
                  // order as in this basis. Check preconditions:
                  UTIL_CHECK(starPtr2->invertFlag == -1);
                  UTIL_CHECK(starId2 = starId + 1);
                  UTIL_CHECK(basisId2 = basisId + 1);
                  UTIL_CHECK(starPtr->waveBz == waveIn);
                  UTIL_CHECK(starPtr2->waveBz == waveIn2);

                  // Copy components for both stars into fields array
                  for (int j = 0; j < nMonomer; ++j) {
                      fields[j][basisId] = temp[j];
                      fields[j][basisId2] = temp2[j];
                  }

               } else
               if (starPtr->invertFlag == -1) {

                  // This is a pair of open stars written in opposite
                  // order from in this basis. Check preconditions:
                  UTIL_CHECK(starPtr2->invertFlag == 1);
                  UTIL_CHECK(starId == starId2 + 1);
                  UTIL_CHECK(basisId == basisId2 + 1);
                  UTIL_CHECK(waveId == starPtr->beginId);

                  // Check that waveIn2 is negation of waveIn
                  IntVec<D> nVec;
                  nVec.negate(waveIn);
                  nVec =
                       shiftToMinimum(nVec, mesh.dimensions(), unitCell);
                  UTIL_CHECK(waveIn2 == nVec);

                  /*
                  * Consider two related stars, C and D, that are listed
                  * in the order (C,D) in the basis used in this code (the
                  * reading program), but that were listed in the opposite
                  * order (D,C) in the program that wrote the file (the
                  * writing program). In the basis of the reading program,
                  * star C has star index starId2, while star D has index
                  * starId = starid2 + 1.
                  *
                  * Let f(r) and f^{*}(r) denote the basis functions used
                  * by the reading program for stars C and D, respectively.
                  * Let u(r) and u^{*}(r) denote the corresponding basis
                  * functions used by the writing program for stars C
                  * and D.  Let exp(i phi) denote the unit magnitude
                  * coefficient (i.e., phasor) within f(r) of the wave
                  * with wave index waveId2, which was the characteristic
                  * wave for star C in the writing program. The
                  * coefficient of this wave within the basis function
                  * u(r) used by the writing program must instead be real
                  * and positive. This implies that
                  * u(r) = exp(-i phi) f(r).
                  *
                  * Focus on the contribution to the field for a specific
                  * monomer type j.  Let a and b denote the desired
                  * coefficients of stars C and D in the reading program,
                  * for which the total contribution of both stars to the
                  * field is:
                  *
                  *  (a - ib) f(r) + (a + ib) f^{*}(r)
                  *
                  * Let A = temp[j] and B = temp2[j] denote the
                  * coefficients read from file in order (A,B).  Noting
                  * that the stars were listed in the order (D,C) in the
                  * basis used by the writing program, the contribution
                  * of both stars must be (A-iB)u^{*}(r)+(A+iB)u(r), or:
                  *
                  *  (A+iB) exp(-i phi)f(r) + (A-iB) exp(i phi) f^{*}(r)
                  *
                  * Comparing coefficients of f^{*}(r), we find that
                  *
                  *       (a + ib) = (A - iB) exp(i phi)
                  *
                  * This equality is implemented below, where the
                  * variable "phasor" is set equal to exp(i phi).
                  */
                  phasor = basis.wave(waveId2).coeff;
                  phasor = phasor/std::abs(phasor);
                  for (int j = 0; j < nMonomer; ++j) {
                      coeff = std::complex<double>(temp[j],-temp2[j]);
                      coeff *= phasor;
                      fields[j][basisId2] = real(coeff);
                      fields[j][basisId ] = imag(coeff);
                  }

                  // Increment count of number of reversed open pairs
                  ++nReversedPair;

               } else {
                  UTIL_THROW("Invalid starInvert value");
               }

            }   // if (wavePtr->invertFlag == 0) ... else ...

         }   // if (waveExists && sizeMatches)

      }   // end while (i < nStarIn)

      if (nReversedPair > 0) {
         Log::file() << "\n";
         Log::file() << nReversedPair << " reversed pairs of open stars"
                     << " detected in readFieldsBasis\n";
      }

   }

   /*
   * Write an array of fields in basis format to an output stream.
   */
   template <int D>
   void writeBasisData(std::ostream &out,
                       DArray<DArray<double> > const & fields,
                       Basis<D> const & basis)
   {
      // Preconditions, set nMonomer and fieldCapacity
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      int fieldCapacity = fields[0].capacity();
      for (int j = 0; j < nMonomer; ++j) {
         UTIL_CHECK(fieldCapacity == fields[j].capacity());
      }
      UTIL_CHECK(basis.isInitialized());

      // Write fields
      int ib = 0; 
      int nStar = basis.nStar();
      for (int i = 0; i < nStar; ++i) {
         if (ib >= fieldCapacity) break;
         if (!basis.star(i).cancel) {
            for (int j = 0; j < nMonomer; ++j) {
               out << Dbl(fields[j][ib], 20, 10);
            }
            out << "   ";
            for (int j = 0; j < D; ++j) {
               out << Int(basis.star(i).waveBz[j], 5);
            }
            out << Int(basis.star(i).size, 5) << std::endl;
            ++ib;
         }
      }

   }

   template <int D, class AT> 
   void convertBasisToKGrid(DArray<double> const & in,
                            AT& out,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions)
   {
      UTIL_CHECK(basis.isInitialized());

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dftDimensions);

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient for star
      std::complex<double> coeff;             // coefficient for wave
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int ib;                                 // basis index
      int iw;                                 // wave index

      // Initialize all dft coponents to zero
      for (rank = 0; rank < dftMesh.size(); ++rank) {
         out[rank][0] = 0.0;
         out[rank][1] = 0.0;
      }

      // Loop over stars, skipping cancelled stars
      is = 0;
      while (is < basis.nStar()) {
         starPtr = &(basis.star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basisId for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            component = std::complex<double>(in[ib], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis.wave(iw);
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Loop over waves in first star
            component = std::complex<double>(in[ib], -in[ib+1]);
            component /= sqrt(2.0);
            starPtr = &(basis.star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis.wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis.star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = std::complex<double>(in[ib], +in[ib+1]);
            component /= sqrt(2.0);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis.wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }

            // Increment is by 2 (two stars were processed)
            is += 2;

         } else {

            UTIL_THROW("Invalid invertFlag value");

         }

      }

   }

   template <int D, class AT>
   void convertKGridToBasis(AT const& in,
                            DArray<double>& out,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions,
                            bool checkSymmetry,
                            double epsilon) 
   {
      UTIL_CHECK(basis.isInitialized());

      if (checkSymmetry) {
         // Check if kgrid has symmetry
         bool symmetric;
         symmetric = hasSymmetry(in, basis, dftDimensions, epsilon, true);
         if (!symmetric) {
            Log::file() << std::endl
               << "WARNING: non-negligible error in conversion to "
               << "symmetry-adapted basis format." << std::endl
               << "   See error values printed above for each "
               << "asymmetric field." << std::endl
               << "   The field that is output by the above operation "
               << "will be a" << std::endl
               << "   symmetrized version of the input field."
               << std::endl << std::endl;
         }
      }

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dftDimensions);

      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // coefficient for star
      int rank;                                // dft grid rank of wave
      int is;                                  // star index
      int ib;                                  // basis index
      int iw;                                  // wave index

      // Initialize all components to zero
      for (is = 0; is < basis.nBasis(); ++is) {
         out[is] = 0.0;
      }

      // Loop over stars
      is = 0;
      while (is < basis.nStar()) {
         starPtr = &(basis.star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basis id for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Choose a wave in the star that is not implicit
            int beginId = starPtr->beginId;
            int endId = starPtr->endId;
            iw = 0;
            bool isImplicit = true;
            while (isImplicit) {
               wavePtr = &basis.wave(beginId + iw);
               if (!wavePtr->implicit) {
                  isImplicit = false;
               } else {
                   UTIL_CHECK(beginId + iw < endId - 1 - iw);
                   wavePtr = &basis.wave(endId - 1 - iw);
                   if (!wavePtr->implicit) {
                      isImplicit = false;
                   }
               }
               ++iw;
            }
            UTIL_CHECK(wavePtr->starId == is);

            // Compute component value
            rank = dftMesh.rank(wavePtr->indicesDft);
            component = std::complex<double>(in[rank][0], in[rank][1]);
            component /= wavePtr->coeff;
            out[ib] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either the first wave of the 1st star or its inverse
            // in the second star
            wavePtr = &basis.wave(starPtr->beginId);
            UTIL_CHECK(wavePtr->starId == is);
            if (wavePtr->implicit) {
               iw = wavePtr->inverseId;
               starPtr = &(basis.star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &basis.wave(iw);
               UTIL_CHECK(!(wavePtr->implicit));
               UTIL_CHECK(wavePtr->starId == is+1);
            }
            rank = dftMesh.rank(wavePtr->indicesDft);
            component = std::complex<double>(in[rank][0], in[rank][1]);
            UTIL_CHECK(std::abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);

            // Compute basis function coefficient values
            if (starPtr->invertFlag == 1) {
               out[ib] = component.real();
               out[ib+1] = -component.imag();
            } else {
               out[ib] = component.real();
               out[ib+1] = component.imag();
            }

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }

   /*
   * Test if an K-grid array has the declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D, class AT>
   bool hasSymmetry(AT const & in, 
                    Basis<D> const& basis,
                    IntVec<D> const& dftDimensions,
                    double epsilon,
                    bool verbose)
   {
      UTIL_CHECK(basis.isInitialized());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> waveCoeff;         // coefficient from wave
      std::complex<double> rootCoeff;         // coefficient from root
      std::complex<double> diff;              // coefficient difference
      int is;                                 // star index
      int iw;                                 // wave index
      int beginId, endId;                     // star begin, end ids
      int rank;                               // dft grid rank of wave

      double cancelledError(0.0);   // max error from cancelled stars
      double uncancelledError(0.0); // max error from uncancelled stars

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dftDimensions);

      // Loop over all stars
      for (is = 0; is < basis.nStar(); ++is) {
         starPtr = &(basis.star(is));

         if (starPtr->cancel) {

            // Check that coefficients are zero for all waves in star
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis.wave(iw);
               if (!wavePtr->implicit) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  waveCoeff 
                     = std::complex<double>(in[rank][0], in[rank][1]);
                  if (std::abs(waveCoeff) > cancelledError) {
                     cancelledError = std::abs(waveCoeff);
                     if ((!verbose) && (cancelledError > epsilon)) {
                        return false;
                     }
                  }
               }
            }

         } else {

            // Check consistency of coeff values from all waves
            bool hasRoot = false;
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis.wave(iw);
               if (!(wavePtr->implicit)) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  waveCoeff 
                     = std::complex<double>(in[rank][0], in[rank][1]);
                  waveCoeff /= wavePtr->coeff;
                  if (hasRoot) {
                     diff = waveCoeff - rootCoeff;
                     if (std::abs(diff) > uncancelledError) {
                        uncancelledError = std::abs(diff);
                        if ((!verbose) && (uncancelledError > epsilon)) {
                           return false;
                        }
                     }
                  } else {
                     rootCoeff = waveCoeff;
                     hasRoot = true;
                  }
               }
            }

         }

      } //  end loop over star index is

      if ((cancelledError < epsilon) && (uncancelledError < epsilon)) {
         return true;
      } else if (verbose) {
         Log::file() << std::endl
                     << "Maximum coefficient of a cancelled star: "
                     << cancelledError << std::endl
                     << "Maximum error of coefficient for uncancelled star: "
                     << uncancelledError << std::endl;
      }
      return false;
   }


} // namespace Prdc
} // namespace Pscf
#endif
