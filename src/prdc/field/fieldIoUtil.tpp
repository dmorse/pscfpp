#ifndef PRDC_FIELD_IO_UTIL_TPP
#define PRDC_FIELD_IO_UTIL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2024, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/fieldHeader.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/complex.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>

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
                            int nMonomer, 
                            IntVec<D> const& dimensions)
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
   * Inspect dimensions of a DArray of fields, each of type FT.
   */
   template <int D, class FT>
   void inspectFields(DArray<FT> const& fields,
                      int & nMonomer,
                      IntVec<D> & dimensions)
   {
      UTIL_CHECK(fields.isAllocated());
      nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      dimensions = fields[0].meshDimensions();
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].meshDimensions() == dimensions);
      }
   }

   /*
   * Check allocation of an array of arrays of type AT, allocate if needed.
   */
   template <class AT>
   void checkAllocateArrays(DArray<AT>& arrays,
                            int nMonomer, 
                            int capacity)
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

   /*
   * Inspect dimensions of an array of arrays of type AT.
   */
   template <class AT>
   void inspectArrays(DArray<AT> const& arrays,
                    int & nMonomer,
                    int & capacity)
   {
      UTIL_CHECK(arrays.isAllocated());
      nMonomer = arrays.capacity();
      UTIL_CHECK(nMonomer > 0);

      capacity = arrays[0].capacity();
      UTIL_CHECK(capacity > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(arrays[i].capacity() == capacity);
      }
   }

   // RGrid File Io Templates

   template <int D>
   void readMeshDimensions(std::istream& in,
                           IntVec<D> const& meshDimensions) 
   {
      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> meshDimensionsIn;
      in >> meshDimensionsIn;
      if (meshDimensionsIn != meshDimensions) {
         Log::file()
           << "Inconsistent mesh dimensions:\n"
           << "meshDimensions (expected)  = " << meshDimensions << "\n"
           << "meshDimensions (from file) = " << meshDimensionsIn << "\n";
         UTIL_THROW("Unexpected mesh dimensions in field file header");
      }
   }

   template <int D>
   void writeMeshDimensions(std::ostream &out,
                            IntVec<D> const& meshDimensions)
   {
      out << "mesh " <<  std::endl
          << "           " << meshDimensions << std::endl;
   }

   template <int D, class ART>
   void readRGridData(std::istream& in,
                      DArray<ART>& fields,
                      int nMonomer,
                      IntVec<D> const& dimensions)
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

   template <int D, class ART>
   void readRGridData(std::istream& in, 
                      ART& field,
                      IntVec<D> const& dimensions)
   {
      MeshIteratorFortran<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> field[rank];
      }
   }

   template <int D, class ART>
   void writeRGridData(std::ostream& out,
                       DArray<ART> const& fields,
                       int nMonomer,
                       IntVec<D> const& dimensions)
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

   template <int D, class ART>
   void writeRGridData(std::ostream& out,
                       ART const& field,
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

   template <int D, class ACT>
   void readKGridData(std::istream& in,
                      DArray<ACT> & fields,
                      int nMonomer,
                      IntVec<D> const& dftDimensions)
   {
      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

      RT x, y;
      MeshIterator<D> iter(dftDimensions);
      int rank, i, j, idum;
      i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == rank);
         for (j = 0; j < nMonomer; ++j) {
            in >> x;
            in >> y;
            assign<CT, RT>(fields[j][rank], x, y);
         }
         ++i;
      }
   }

   template <int D, class ACT>
   void readKGridData(std::istream& in, 
                      ACT& field,
                      IntVec<D> const& dftDimensions)
   {
      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

      RT x, y;
      MeshIterator<D> iter(dftDimensions);
      int rank, idum;
      int i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == rank);
         //in >> field[rank][0];
         //in >> field[rank][1];
         in >> x;
         in >> y;
         assign<CT, RT>(field[rank], x, y);
         ++i;
      }
   }

   template <int D, class ACT>
   void writeKGridData(std::ostream& out,
                       DArray<ACT> const& fields,
                       int nMonomer,
                       IntVec<D> const& dftDimensions)
   {
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());

      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

      RT x, y;
      MeshIterator<D> iter(dftDimensions);
      int rank;
      int i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         UTIL_CHECK(i == rank);
         out << Int(rank, 5);
         for (int j = 0; j < nMonomer; ++j) {
            x = real<CT, RT>(fields[j][rank]);
            y = imag<CT, RT>(fields[j][rank]);
            out << "  "
                << Dbl(x, 21, 13)
                << Dbl(y, 21, 13);
         }
         out << std::endl;
         ++i;
      }

   }

   template <int D, class ACT>
   void writeKGridData(std::ostream& out, 
                       ACT const& field,
                       IntVec<D> const& dftDimensions)
   {
      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

      RT x, y;
      MeshIterator<D> iter(dftDimensions);
      int rank, i;
      i = 0;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         UTIL_CHECK(i == rank);
         x = real<CT, RT>(field[rank]);
         y = imag<CT, RT>(field[rank]);
         out << Int(rank, 5);
         out << "  "
             << Dbl(x, 21, 13)
             << Dbl(y, 21, 13);
         out << std::endl;
         ++i;
      }
   }

   // Functions for files in symmetry-adapted basis format

   /*
   * Read the number of basis functions from a field file header.
   */
   int readNBasis(std::istream& in)
   {
   
      // Read the label, which can be N_star or N_basis
      std::string label;
      in >> label;
      if (label != "N_star" && label != "N_basis") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected N_basis or N_star, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
 
      // Read the value of nBasis
      int nBasis;
      in >> nBasis;
      UTIL_CHECK(nBasis > 0);

      return nBasis;
   }

   /*
   * Write the number of basis functions to a field file header.
   */
   void writeNBasis(std::ostream& out, int nBasis)
   {
      out << "N_basis      " << std::endl
          << "             " << nBasis << std::endl;
   }

   /*
   * Read the data section for an array of fields in basis format.
   */
   template <int D>
   void readBasisData(std::istream& in,
                      DArray< DArray<double> >& fields,
                      UnitCell<D> const& unitCell,
                      Mesh<D> const& mesh,
                      Basis<D> const& basis,
                      int nBasisIn)
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

      // Reset nBasis = min(nBasisIn, fieldCapacity)
      int nBasis = nBasisIn;
      if (fieldCapacity < nBasisIn) {
         nBasis = fieldCapacity;
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
      while (i < nBasis) {

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

      }   // end while (i < nBasis)

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
                       DArray< DArray<double> > const & fields,
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

   template <int D, class ACT> 
   void convertBasisToKGrid(DArray<double> const & in,
                            ACT& out,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions)
   {
      UTIL_CHECK(basis.isInitialized());

      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

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
         //out[rank][0] = 0.0;
         //out[rank][1] = 0.0;
         assign<CT, RT>(out[rank], 0.0, 0.0);
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
                  //out[rank][0] = coeff.real();
                  //out[rank][1] = coeff.imag();
                  assign<CT, RT>(out[rank], coeff);
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
                  //out[rank][0] = coeff.real();
                  //out[rank][1] = coeff.imag();
                  assign<CT, RT>(out[rank], coeff);
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
                  //out[rank][0] = coeff.real();
                  //out[rank][1] = coeff.imag();
                  assign<CT, RT>(out[rank], coeff);
               }
            }

            // Increment is by 2 (two stars were processed)
            is += 2;

         } else {

            UTIL_THROW("Invalid invertFlag value");

         }

      }

   }

   template <int D, class ACT>
   void convertKGridToBasis(ACT const& in,
                            DArray<double>& out,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions,
                            bool checkSymmetry,
                            double epsilon) 
   {
      UTIL_CHECK(basis.isInitialized());

      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

      // Check if input field in k-grid format has specified symmetry
      if (checkSymmetry) {
         bool symmetric;
         symmetric = hasSymmetry(in, basis, dftDimensions, epsilon, true);
         if (!symmetric) {
            Log::file() << std::endl
               << "WARNING: Conversion of asymmetric field to"
               << "symmetry-adapted basis in Prdc::convertKgridToBasis." 
               << std::endl << std::endl;
         }
      }

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dftDimensions);

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient for star
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int ib;                                 // basis index
      int iw;                                 // wave index

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
            //component = std::complex<double>(in[rank][0], in[rank][1]);
            assign<CT,RT>(component, in[rank]);
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
            //component = std::complex<double>(in[rank][0], in[rank][1]);
            assign<CT, RT>(component, in[rank]); 
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
   template <int D, class ACT>
   bool hasSymmetry(ACT const & in, 
                    Basis<D> const& basis,
                    IntVec<D> const& dftDimensions,
                    double epsilon,
                    bool verbose)
   {
      UTIL_CHECK(basis.isInitialized());

      typedef typename ACT::Complex CT;
      typedef typename ACT::Real    RT;

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
                  //waveCoeff 
                  //   = std::complex<double>(in[rank][0], in[rank][1]);
                  assign<CT, RT>(waveCoeff, in[rank]);
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
                  //waveCoeff 
                  //   = std::complex<double>(in[rank][0], in[rank][1]);
                  assign<CT, RT>(waveCoeff, in[rank]);
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
         Log::file() 
              << std::endl
              << "Maximum coefficient of a cancelled star: "
              << cancelledError << std::endl
              << "Maximum error of coefficient for an uncancelled star: "
              << uncancelledError << std::endl;
      }
      return false;
   }

   // Field file manipulations

   template <int D, class ART>
   void replicateUnitCell(std::ostream &out,
                          DArray<ART> const & fields,
                          IntVec<D> const & meshDimensions,
                          UnitCell<D> const & unitCell,
                          IntVec<D> const & replicas)
   {
      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Compute properties of mesh for replicated fields
      IntVec<D> repDimensions;  // Dimensions
      int repSize = 1;          // Total number of grid points
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(replicas[i] > 0);
         UTIL_CHECK(meshDimensions[i] > 0);
         repDimensions[i] = replicas[i] * meshDimensions[i];
         repSize *= repDimensions[i];
      }

      // Allocate arrays for replicated fields
      DArray< DArray<double> > repFields;
      repFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         repFields[i].allocate(repSize);
      }

      // Replicate data
      Mesh<D> mesh(meshDimensions);
      Mesh<D> cellMesh(replicas);
      Mesh<D> repMesh(repDimensions);
      MeshIterator<D> iter(meshDimensions);
      MeshIterator<D> cellIter(replicas);
      IntVec<D> position;
      IntVec<D> cellPosition;
      IntVec<D> repPosition;
      double value;
      int repRank;
      for (int i = 0; i < nMonomer; ++i) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            position = iter.position();
            value = fields[i][iter.rank()];
            for (cellIter.begin(); !cellIter.atEnd(); ++cellIter) {
               cellPosition = cellIter.position(); 
               for (int j=0; j < D; ++j) {
                  repPosition[j] = position[j] 
                                 + meshDimensions[j]*cellPosition[j];
               }
               repRank = repMesh.rank(repPosition);
               repFields[i][repRank] = value;
            }
         }
      }

      // Set up UnitCell for replicated fields
      UnitCell<D> cell;
      FSArray<double, 6> parameters;
      int nParameter = unitCell.nParameter();
      for (int i = 0; i < nParameter; i++) {
         parameters[i]=  replicas[i]* unitCell.parameter(i);
      }
      cell.set(unitCell.lattice(), parameters);

      // Write header
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
      writeMeshDimensions(out, repDimensions);

      // Write field data
      writeRGridData(out, repFields, nMonomer, repDimensions);
   }

   template <int D, class ART>
   void
   expandRGridDimension(std::ostream &out,
                        DArray<ART> const & fields,
                        IntVec<D> const & meshDimensions,
                        UnitCell<D> const & unitCell,
                        int d,
                        DArray<int> newGridDimensions)
   {
      UTIL_THROW("Unimplemented base template"); 
   }

   template <class ART>
   void
   expandRGridDimension(std::ostream &out,
                        DArray<ART> const & fields,
                        IntVec<1> const & meshDimensions,
                        UnitCell<1> const & unitCell,
                        int d,
                        DArray<int> newGridDimensions)

   {
      // Check validity of expanded dimension d and newGridDimensions
      UTIL_CHECK(d > 1);
      UTIL_CHECK(d <= 3);
      UTIL_CHECK(newGridDimensions.capacity() == (d - 1));


      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
 
      // Set up necessary objects
      FSArray<double, 6> cellParameters;
      cellParameters.append(unitCell.parameter(0));
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
 
      if (d == 2) {

         // 1D expanded to 2D

         // Set dimensions
         IntVec<2> dimensions;
         dimensions[0] = meshDimensions[0];
         dimensions[1] = newGridDimensions[0];

         // Assign unit cell
         UnitCell<2> cell;
         if (dimensions[0] == dimensions[1]) {
            cell.set(UnitCell<2>::Square, cellParameters);
         } else {
            cellParameters.append((double)dimensions[1]/dimensions[0]
                                                  * cellParameters[0]);
            cell.set(UnitCell<2>::Rectangular, cellParameters);
         }

         // Write Header
         Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
         out << "mesh " <<  std::endl
             << "           " << dimensions << std::endl;

      } else if (d == 3) {

         // 1D expanded to 3D

         // Set dimensions
         IntVec<3> dimensions;
         dimensions[0] = meshDimensions[0];
         dimensions[1] = newGridDimensions[0];
         dimensions[2] = newGridDimensions[1];

         // Assign unit cell
         UnitCell<3> cell;
         if (dimensions[2] == dimensions[1]) {
            if (dimensions[1] == dimensions[0]) {
               cell.set(UnitCell<3>::Cubic, cellParameters);
            } else {
               cellParameters.append((double)dimensions[1]/dimensions[0]
                                                     * cellParameters[0]);
               cellParameters.append((double)dimensions[2]/dimensions[0]
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Orthorhombic, cellParameters);
            }
         } else {
            if (dimensions[1] == dimensions[0]) {
               cellParameters.append((double)dimensions[2]/dimensions[0]
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Tetragonal, cellParameters);
            } else {
               cellParameters.append((double)dimensions[1]/dimensions[0]
                                                     * cellParameters[0]);
               cellParameters.append((double)dimensions[2]/dimensions[0]
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Orthorhombic, cellParameters);
            }
         }

         // Write Header
         Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
         out << "mesh " <<  std::endl
             << "           " << dimensions << std::endl;

      } else {

         UTIL_THROW("Invalid d value");

      }

      // Write expanded fields
      int nReplica = newGridDimensions[0];
      if (d == 3) {
         nReplica *= newGridDimensions[1];
      }
      MeshIteratorFortran<1> iter(meshDimensions);
      int rank;
      for (int i = 0; i < nReplica; ++i) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(fields[j][rank], 21, 13);
            }
            out << std::endl;
         }
      }

   }

   template <class ART>
   void expandRGridDimension(std::ostream &out,
                             DArray<ART> const & fields,
                             IntVec<2> const & meshDimensions,
                             UnitCell<2> const & unitCell,
                             int d,
                             DArray<int> newGridDimensions)
   {
      // 2D expanded to 3D

      // Check validity of expanded dimension d and newGridDimensions
      UTIL_CHECK(d == 3);
      UTIL_CHECK(newGridDimensions.capacity() == 1);

      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Set dimensions of mesh for replicated fields
      IntVec<3> dimensions;
      dimensions[0] = meshDimensions[0];
      dimensions[1] = meshDimensions[1];
      dimensions[2] = newGridDimensions[0];

      // Set unit cell for replicated fields
      UnitCell<3> cell;
      FSArray<double, 6> cellParameters;
      cellParameters.append(unitCell.parameter(0));
      if (unitCell.lattice() == UnitCell<2>::Square) {
         if (newGridDimensions[0] == meshDimensions[0]){
            cell.set(UnitCell<3>::Cubic, cellParameters);
         } else {
            cellParameters.append((double)dimensions[2]/dimensions[0]
                                                  * cellParameters[0]);
            cell.set(UnitCell<3>::Tetragonal, cellParameters);
         }
      } else if (unitCell.lattice() == UnitCell<2>::Rectangular) {
         cellParameters.append(unitCell.parameter(1));
         cellParameters.append((double)dimensions[2]/dimensions[0]
                                               * cellParameters[0]);
         cell.set(UnitCell<3>::Orthorhombic, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Hexagonal){
         cellParameters.append((double)dimensions[2]/dimensions[0]
                                               * cellParameters[0]);
         cell.set(UnitCell<3>::Hexagonal, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Rhombic) {
         cellParameters.append(unitCell.parameter(0));
         cellParameters.append((double)dimensions[2]/dimensions[0]
                                               * cellParameters[0]);
         cellParameters.append(Constants::Pi / 2);
         cellParameters.append(0.0);
         cellParameters.append(unitCell.parameter(1));
         cell.set(UnitCell<3>::Triclinic, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Oblique) {
         cellParameters.append(unitCell.parameter(1));
         cellParameters.append((double)dimensions[2]/dimensions[0]
                                               * cellParameters[0]);
         cellParameters.append(Constants::Pi / 2);
         cellParameters.append(0.0);
         cellParameters.append(unitCell.parameter(2));
         cell.set(UnitCell<3>::Triclinic, cellParameters);
      } else {
         UTIL_THROW("Unrecognized 2D lattice system.");
      }

      // Write Header
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
      out << "mesh " <<  std::endl
          << "           " << dimensions << std::endl;

      // Write expanded fields
      int nReplica = newGridDimensions[0];
      MeshIteratorFortran<2> iter(meshDimensions);
      int rank;
      for (int i = 0; i < nReplica; ++i) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(fields[j][rank], 21, 13);
            }
            out << std::endl;
         }
      }

   }

   template <class ART>
   void expandRGridDimension(std::ostream &out,
                             DArray<ART> const & fields,
                             IntVec<3> const & meshDimensions,
                             UnitCell<3> const & unitCell,
                             int d,
                             DArray<int> newGridDimensions)

   {  UTIL_THROW("expandRGridDimension is invalid when D = 3."); }

} // namespace Prdc
} // namespace Pscf
#endif
