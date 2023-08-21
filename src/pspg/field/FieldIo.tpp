#ifndef PSPG_FIELD_IO_TPP
#define PSPG_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/misc/Log.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>
#include <cmath>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : meshPtr_(0),
      fftPtr_(0),
      latticePtr_(0),
      groupNamePtr_(0),
      groupPtr_(0),
      basisPtr_(0),
      fileMasterPtr_(0)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldIo<D>::~FieldIo()
   {}

   /*
   * Get and store addresses of associated objects.
   */
   template <int D>
   void FieldIo<D>::associate(Mesh<D> const & mesh,
                              FFT<D> const & fft,
                              typename UnitCell<D>::LatticeSystem & lattice,
                              std::string& groupName,
                              SpaceGroup<D>& group,
                              Basis<D>& basis,
                              FileMaster const & fileMaster)
   {
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      latticePtr_ = &lattice;
      groupNamePtr_ = &groupName;
      groupPtr_ = &group;
      basisPtr_ = &basis;
      fileMasterPtr_ = &fileMaster;
   }

   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in,
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      // Read generic part of field file header.
      // Set nMonomer to number of field components in file
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);
      UTIL_CHECK(nMonomer > 0);

      // Read number of stars from field file into StarIn
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
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      // Check dimensions of fields array
      int nMonomerFields = fields.capacity();
      UTIL_CHECK(nMonomer == nMonomerFields);
      int i, j;
      int fieldCapacity = fields[0].capacity();
      for (i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].capacity() == fieldCapacity);
      }

      // Define nStar and nBasis
      int nStar = basis().nStar();
      int nBasis = basis().nBasis();

      // Initialize all elements of field components to zero
      for (j = 0; j < nMonomer; ++j) {
         UTIL_CHECK(fields[j].capacity() == nBasis);
         for (i = 0; i < nBasis; ++i) {
            fields[j][i] = 0.0;
         }
      }

      // Allocate buffer arrays used to read components
      DArray<double> temp, temp2;
      temp.allocate(nMonomer);
      temp2.allocate(nMonomer);

      // Variables needed in loop over stars
      typename Basis<D>::Star const * starPtr;
      typename Basis<D>::Star const * starPtr2;
      IntVec<D> waveIn, waveIn2;
      int starId, starId2;
      int basisId, basisId2;
      int waveId, waveId2;

      std::complex<double> coeff, phasor;
      IntVec<D> waveBz, waveDft;
      int nWaveVector;
      int nReversedPair = 0;
      bool waveExists;

      // Loop over stars in input file to read field components
      i = 0;
      while (i < nStarIn) {

         // Read next line of data
         for (int j = 0; j < nMonomer; ++j) {
            in >> temp[j];               // field components
         }
         in >> waveIn;                   // wave of star
         in >> nWaveVector;              // # of waves in star

         // Check if waveIn is in first Brillouin zone (FBZ) for the mesh.
         waveBz = shiftToMinimum(waveIn, mesh().dimensions(), unitCell);
         waveExists = (waveIn == waveBz);

         if (!waveExists) {

            //  If wave is not in FBZ, ignore and continue
            ++i;

         } else {

            // If wave is in FBZ, process the line

            // Find the star containing waveIn
            waveDft = waveIn;
            mesh().shift(waveDft);
            waveId = basis().waveId(waveDft);
            starId = basis().wave(waveId).starId;
            starPtr = &basis().star(starId);
            UTIL_CHECK(!(starPtr->cancel));
            // basisId = starId;
            basisId = starPtr->basisId;

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
                     <<  "waveBz of star = " << starPtr->waveBz  << "\n";
               }
               ++i;  // increment input line counter i

            } else {       // If starPtr->invertFlag != 0


               // Read the next line
               for (int j = 0; j < nMonomer; ++j) {
                  in >> temp2[j];               // components of field
               }
               in >> waveIn2;                   // wave of star
               in >> nWaveVector;               // # of wavevectors in star

               // Check that waveIn2 is also in the 1st BZ
               waveBz =
                   shiftToMinimum(waveIn2, mesh().dimensions(), unitCell);
               UTIL_CHECK(waveIn2 == waveBz);

               // Identify the star containing waveIn2
               waveDft = waveIn2;
               mesh().shift(waveDft);
               waveId2 = basis().waveId(waveDft);
               starId2 = basis().wave(waveId2).starId;
               starPtr2 = &basis().star(starId2);
               UTIL_CHECK(!(starPtr2->cancel));
               //basisId2 = starId2;
               basisId2 = starPtr2->basisId;

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
                       shiftToMinimum(nVec, mesh().dimensions(), unitCell);
                  UTIL_CHECK(waveIn2 == nVec);

                  /*
                  * Consider two related stars, C and D, that are listed in
                  * the order (C,D) in the basis used in this code (the
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
                  phasor = basis().wave(waveId2).coeff;
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

               // Increment counter by 2 because two lines were read
               i = i + 2;

            }   // if (wavePtr->invertFlag == 0) ... else ...
         }   // if (!waveExists) ... else ...
      }   // end while (i < nStarIn)

      if (nReversedPair > 0) {
         Log::file() << "\n";
         Log::file() << nReversedPair << " reversed pairs of open stars"
                     << " detected in FieldIo::readFieldsBasis\n";
      }

   }


   template <int D>
   void FieldIo<D>::readFieldsBasis(std::string filename,
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
       std::ifstream file;
       fileMaster().openInputFile(filename, file);
       readFieldsBasis(file, fields, unitCell);
       file.close();
   }

   template <int D>
   void 
   FieldIo<D>::writeFieldsBasis(std::ostream &out,
                                DArray< DArray<double> > const&  fields,
                                UnitCell<D> const & unitCell)
   const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(basis().isInitialized());

      // Write header
      writeFieldHeader(out, nMonomer, unitCell);
      out << "N_basis      " << std::endl
          << "             " << basis().nBasis() << std::endl;

      // Define nStar and nBasis
      int nStar = basis().nStar();
      //int nBasis = nStar;
      int nBasis = basis().nBasis();

      // Write fields
      int ib = 0;
      for (int i = 0; i < nStar; ++i) {
         if (!basis().star(i).cancel) {
            UTIL_CHECK(ib == basis().star(i).basisId);
            for (int j = 0; j < nMonomer; ++j) {
               out << Dbl(fields[j][ib], 20, 10);
            }
            out << "   ";
            for (int j = 0; j < D; ++j) {
               out << Int(basis().star(i).waveBz[j], 5);
            }
            out << Int(basis().star(i).size, 5) << std::endl;
            ++ib;
         }
      }

   }

   template <int D>
   void 
   FieldIo<D>::writeFieldsBasis(std::string filename,
                                DArray< DArray<double> > const& fields,
                                UnitCell<D> const & unitCell)
   const
   {
       std::ofstream file;
       fileMaster().openOutputFile(filename, file);
       writeFieldsBasis(file, fields, unitCell);
       file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::istream &in,
                                    DArray<RDField<D> >& fields,
                                    UnitCell<D>& unitCell) const
   {

      // Read generic part of field header
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());

      // Read grid dimensions
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
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      DArray<cudaReal*> temp_out;
      temp_out.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaReal[mesh().size()];
      }

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         for(int k = 0; k < nMonomer; ++k) {
            in >> std::setprecision(15)>> temp_out[k][rank];
         }
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         }
      }

      for (int i = 0; i < nMonomer; i++) {
         cudaMemcpy(fields[i].cDField(), temp_out[i],
            mesh().size() * sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename,
                                    DArray< RDField<D> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields, unitCell);
      file.close();
   }
   
   template <int D>
   void FieldIo<D>::readFieldsRGridData(std::istream &in,
                                        DArray< RDField<D> >& fields,
                                        int nMonomer) const
   {
      // Read grid dimensions
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
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      DArray<cudaReal*> temp_out;
      temp_out.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaReal[mesh().size()];
      }

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         for(int k = 0; k < nMonomer; ++k) {
            in >> std::setprecision(15)>> temp_out[k][rank];
         }
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         }
      }

      for (int i = 0; i < nMonomer; i++) {
         cudaMemcpy(fields[i].cDField(), temp_out[i],
            mesh().size() * sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::ostream &out,
                                     DArray<RDField<D> > const& fields,
                                     UnitCell<D> const & unitCell,
                                     bool writeHeader)
   const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell);
      }
      out << "mesh " <<  std::endl
          << "           " << mesh().dimensions() << std::endl;

      DArray<cudaReal*> temp_out;
      temp_out.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaReal[mesh().size()];
         cudaMemcpy(temp_out[i], fields[i].cDField(),
                    mesh().size() * sizeof(cudaReal), cudaMemcpyDeviceToHost);
      }

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         for(int k = 0; k < nMonomer; ++k) {
            out << "  " << Dbl(temp_out[k][rank], 18, 15);
         }
         out<<'\n';
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         }
      }

      for(int i = 0; i < nMonomer; ++i) {
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }

   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::string filename,
                                     DArray< RDField<D> > const & fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::istream &in,
                                   RDField<D> &field,
                                   UnitCell<D>& unitCell)
   const
   {
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);
      UTIL_CHECK(nMonomer == 1);

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
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      cudaReal* temp_out = new cudaReal[mesh().size()];

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         in >> std::setprecision(15) >> temp_out[rank];

         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         }
      }

      cudaMemcpy(field.cDField(), temp_out,
            mesh().size() * sizeof(cudaReal), cudaMemcpyHostToDevice);
      delete[] temp_out;
      temp_out = nullptr;

   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::string filename,
                                   RDField<D> &field,
                                   UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldRGrid(file, field, unitCell);
      file.close();
   }
   
   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::ostream &out,
                                    RDField<D> const & field,
                                    UnitCell<D> const & unitCell,
                                    bool writeHeader)
   const
   {

      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell);
         out << "mesh " <<  std::endl
             << "           " << mesh().dimensions() << std::endl;
      }

      cudaReal* temp_out = new cudaReal[mesh().size()];;
      cudaMemcpy(temp_out, field.cDField(),
                  mesh().size() * sizeof(cudaReal), cudaMemcpyDeviceToHost);

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         out << "  " << Dbl(temp_out[rank], 18, 15);
         out<<'\n';
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         }
      }

      delete[] temp_out;
      temp_out = nullptr;
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::string filename,
                                    RDField<D> const & field,
                                    UnitCell<D> const & unitCell)
   const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::istream &in,
                                    DArray< RDFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {

      // Read header
      int nMonomer;
      readFieldHeader(in, nMonomer, unitCell);
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nMonomer == fields.capacity());
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
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      DArray<cudaComplex*> temp_out;
      temp_out.allocate(nMonomer);

      int kSize = 1;
      for (int i = 0; i < D; i++) {
         if (i == D - 1) {
            kSize *= (mesh().dimension(i) / 2 + 1);
         }
         else {
            kSize *= mesh().dimension(i);
         }
      }

      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaComplex[kSize];
      }

      // Read Fields;
      int idum;
      MeshIterator<D> itr(mesh().dimensions());
      for (int i = 0; i < kSize; ++i) {
         in >> idum;
         for (int j = 0; j < nMonomer; ++j) {
            in >> temp_out [j][i].x;
            in >> temp_out [j][i].y;
         }
      }

      for(int i = 0; i < nMonomer; ++i) {
         cudaMemcpy(fields[i].cDField(), temp_out[i],
            kSize * sizeof(cudaComplex), cudaMemcpyHostToDevice);
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::string filename,
                                    DArray< RDFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::ostream &out,
                                     DArray<RDFieldDft<D> > const& fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Write header
      writeFieldHeader(out, nMonomer, unitCell);
      out << "mesh " << std::endl
          << "               " << mesh().dimensions() << std::endl;

      DArray<cudaComplex*> temp_out;
      int kSize = 1;
      for (int i = 0; i < D; i++) {
         if (i == D - 1) {
            kSize *= (mesh().dimension(i) / 2 + 1);
         }
         else {
            kSize *= mesh().dimension(i);
         }
      }
      temp_out.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaComplex[kSize];
         cudaMemcpy(temp_out[i], fields[i].cDField(),
            kSize * sizeof(cudaComplex), cudaMemcpyDeviceToHost);
      }

      // Write fields
      MeshIterator<D> itr(mesh().dimensions());
      for (int i = 0; i < kSize; i++) {
         out << Int(i, 5);
         for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(temp_out[j][i].x, 18, 11)
              << Dbl(temp_out[j][i].y, 18, 11);
         }
         out << std::endl;
      }

      for(int i = 0; i < nMonomer; ++i) {
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::string filename,
                                     DArray< RDFieldDft<D> > const & fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   #if 0
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in,
                                    int& nMonomer,
                                    UnitCell<D>& unitCell) const
   {
      // Preconditions
      UTIL_CHECK(latticePtr_);
      UTIL_CHECK(groupNamePtr_);

      int ver1, ver2;
      std::string groupNameIn;
      Pscf::readFieldHeader(in, ver1, ver2, unitCell, groupNameIn, nMonomer);
      // Note:: Function definition in pscf/crystal/UnitCell.tpp

      // Check data from header
      UTIL_CHECK(ver1 == 1);
      UTIL_CHECK(ver2 == 0);
      UTIL_CHECK(unitCell.isInitialized());
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.nParameter() > 0);

      UTIL_CHECK(nMonomer > 0);
      if (groupNameIn != groupName()) {
         Log::file() << std::endl
             << "Warning - "
             << "Mismatched group names in FieldIo::readFieldHeader \n"
             << "  FieldIo::groupName   :" << groupName() << "\n"
             << "  Field file groupName :" << groupNameIn << "\n";
      }

   }
   #endif

   /*
   * Read common part of field header.
   *
   * Return the number of monomers in the file.
   * Set the unitCell.
   * If necessary, make the Basis.
   */
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in, 
                                    int& nMonomer,
                                    UnitCell<D>& unitCell) 
   const
   {
      // Preconditions
      UTIL_CHECK(latticePtr_);
      UTIL_CHECK(groupNamePtr_);
      if (unitCell.lattice() == UnitCell<D>::Null) {
         UTIL_CHECK(unitCell.nParameter() == 0);
      } else {
         UTIL_CHECK(unitCell.nParameter() > 0);
         if (lattice() != UnitCell<D>::Null) {
            UTIL_CHECK(unitCell.lattice() == lattice());
         }
      }

      // Read field header to set unitCell, groupNameIn, nMonomer
      int ver1, ver2;
      std::string groupNameIn;

      Pscf::readFieldHeader(in, ver1, ver2, unitCell, 
                            groupNameIn, nMonomer);
      // Note: Function definition in pscf/crystal/UnitCell.tpp

      // Checks of data from header
      UTIL_CHECK(ver1 == 1);
      UTIL_CHECK(ver2 == 0);
      UTIL_CHECK(unitCell.isInitialized());
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.nParameter() > 0);

      // Validate or initialize lattice type
      if (lattice() == UnitCell<D>::Null) {
         lattice() = unitCell.lattice();
      } else {
         if (lattice() != unitCell.lattice()) {
            Log::file() << std::endl 
               << "Error - "
               << "Mismatched lattice types, FieldIo::readFieldHeader:\n" 
               << "  FieldIo::lattice  :" << lattice() << "\n"
               << "  Unit cell lattice :" << unitCell.lattice() 
               << "\n";
            UTIL_THROW("Mismatched lattice types");
         }
      }

      // Validate or initialize group name
      if (groupName() == "") {
         groupName() = groupNameIn;
      } else {
         if (groupNameIn != groupName()) {
            Log::file() << std::endl 
               << "Error - "
               << "Mismatched group names in FieldIo::readFieldHeader:\n" 
               << "  FieldIo::groupName :" << groupName() << "\n"
               << "  Field file header  :" << groupNameIn << "\n";
            UTIL_THROW("Mismatched group names");
         }
      }

      // Check group, read from file if necessary
      UTIL_CHECK(groupPtr_);
      if (group().size() == 1) {
         if (groupName() != "I") {
            readGroup(groupName(), group());
         }
      }

      // Check basis, construct if not initialized
      UTIL_CHECK(basisPtr_);
      if (!basis().isInitialized()) {
         basisPtr_->makeBasis(mesh(), unitCell, group());
      }
      UTIL_CHECK(basis().isInitialized());
      
   }


   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out,
                                     int nMonomer,
                                     UnitCell<D> const & unitCell) const
   {
      out << "format  1   0" <<  std::endl;
      out << "dim" <<  std::endl
          << "          " << D << std::endl;
      writeUnitCellHeader(out, unitCell);
      out << "group_name" << std::endl
          << "          " << groupName() <<  std::endl;
      out << "N_monomer"  << std::endl
          << "          " << nMonomer << std::endl;
   }

   template <int D>
   void 
   FieldIo<D>::convertBasisToKGrid(DArray<double> const& components,
                                   RDFieldDft<D>& dft) const
   {
      UTIL_CHECK(mesh().dimensions() == dft.meshDimensions());
      UTIL_CHECK(basis().isInitialized());
      const int nBasis = basis().nBasis();
      const int nStar = basis().nStar();
      UTIL_CHECK(nBasis > 0);
      UTIL_CHECK(nStar >= nBasis);

      // Compute DFT array size
      int kSize = 1;
      for (int i = 0; i < D; i++) {
         if (i == D - 1) {
            kSize *= (mesh().dimension(i) / 2 + 1);
         }
         else {
            kSize *= mesh().dimension(i);
         }
      }
      kSize = dft.capacity();
      
      // Allocate dft_out
      cudaComplex* dft_out;
      dft_out = new cudaComplex[kSize];


      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient for star
      std::complex<double> coeff;             // coefficient for wave
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int ib;                                 // basis index
      int iw;                                 // wave index

      // Initialize all dft_out components to zero
      for (rank = 0; rank < dftMesh.size(); ++rank) {
         dft_out[rank].x = 0.0;
         dft_out[rank].y = 0.0;
      }

      // Loop over stars, skipping cancelled stars
      is = 0;
      while (is < basis().nStar()) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basis id for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            UTIL_ASSERT(!isnan(components[ib]));
            UTIL_ASSERT(!isinf(components[ib]));
            component = std::complex<double>(components[ib], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  UTIL_ASSERT(!isnan(coeff.real()));
                  UTIL_ASSERT(!isnan(coeff.imag()));
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = coeff.real();
                  dft_out[rank].y = coeff.imag();
               }
            }
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            UTIL_ASSERT(!isnan(components[ib]));
            UTIL_ASSERT(!isnan(components[ib+1]));
            UTIL_ASSERT(!isinf(components[ib]));
            UTIL_ASSERT(!isinf(components[ib+1]));

            // Loop over waves in first star
            component = std::complex<double>(components[ib],
                                             -components[ib+1]);
            component /= sqrt(2.0);
            starPtr = &(basis().star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  UTIL_ASSERT(!isnan(coeff.real()));
                  UTIL_ASSERT(!isnan(coeff.imag()));
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = (cudaReal)coeff.real();
                  dft_out[rank].y = (cudaReal)coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis().star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = std::complex<double>(components[ib],
                                             +components[ib+1]);
            component /= sqrt(2.0);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  UTIL_ASSERT(!isnan(coeff.real()));
                  UTIL_ASSERT(!isnan(coeff.imag()));
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = (cudaReal)coeff.real();
                  dft_out[rank].y = (cudaReal)coeff.imag();
               }
            }

            // Increment star counter by 2 (two stars were processed)
            is += 2;

         } else {

            UTIL_THROW("Invalid invertFlag value");

         }
      }

      // Copy dft_out (host) to dft (device)
      cudaMemcpy(dft.cDField(), dft_out,
                 kSize * sizeof(cudaComplex), cudaMemcpyHostToDevice);

      delete[] dft_out;

   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(RDFieldDft<D> const& dft,
                                        DArray<double>& components) 
   const
   {
      UTIL_CHECK(basis().isInitialized());
      const int nStar = basis().nStar();       // Number of stars
      const int nBasis = basis().nBasis();     // Number of basis functions

      // Allocate dft_in
      int kSize = 1;
      for (int i = 0; i < D; i++) {
         if (i == D - 1) {
            kSize *= (mesh().dimension(i) / 2 + 1);
         }
         else {
            kSize *= mesh().dimension(i);
         }
      }
      cudaComplex* dft_in;
      dft_in = new cudaComplex[kSize];

      // Copy dft (device) to dft_in (host)
      cudaMemcpy(dft_in, dft.cDField(),
             kSize * sizeof(cudaComplex), cudaMemcpyDeviceToHost);

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      // Declare variables needed in loop over stars
      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // component for star
      std::complex<double> coefficient;        // coefficient in basis
      int rank;                                // dft grid rank of wave
      int is;                                  // star index
      int ib;                                  // basis index
      int iw;                                  // wave id, within star
      bool isImplicit;

      // Initialize all elements of components to zero
      for (is = 0; is < nBasis; ++is) {
         components[is] = 0.0;
      }

      // Loop over all stars 
      is = 0;
      while (is < nStar) {
         starPtr = &(basis().star(is));

         // Skip if star is cancelled
         if (starPtr->cancel) {
            ++is;
            continue;
         }

         // Set basis id for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Choose a characteristic wave that is not implicit.
            // Start with the first, alternately searching from
            // the beginning and end of star.
            int beginId = starPtr->beginId;
            int endId = starPtr->endId;
            iw = 0;
            isImplicit = true;
            while (isImplicit) {
                wavePtr = &basis().wave(beginId + iw);
                if (!wavePtr->implicit) {
                   isImplicit = false;
                } else {
                   UTIL_CHECK(beginId + iw < endId - 1 - iw);
                   wavePtr = &basis().wave(endId - 1 - iw);
                   if (!wavePtr->implicit) {
                      isImplicit = false;
                   }
                }
                ++iw;
            }
            UTIL_CHECK(wavePtr->starId == is);

            // Compute component value 
            rank = dftMesh.rank(wavePtr->indicesDft);
            component =
                   std::complex<double>(dft_in[rank].x, dft_in[rank].y);
            bool componentError = false;
            #if UTIL_DEBUG
            if (std::isnan(component.real())) componentError = true;
            if (std::isnan(component.imag())) componentError = true;
            #endif
            coefficient = wavePtr->coeff;
            #if UTIL_DEBUG
            if (std::isnan(coefficient.real())) componentError = true;
            if (std::isnan(coefficient.imag())) componentError = true;
            if (std::isnormal(coefficient))     componentError = true;
            #endif
            component /= coefficient;
            #if UTIL_DEBUG
            if (std::isnan(component.real())) componentError = true;
            if (std::isnan(component.imag())) componentError = true;
            if (std::isinf(component.real())) componentError = true;
            if (std::isinf(component.imag())) componentError = true;
            #endif

            // Verify that imaginary component is very small
            #ifdef SINGLE_PRECISION
            if (abs(component.imag()) > 1.0E-4) componentError = true;
            #else
            if (abs(component.imag()) > 1.0E-8) componentError = true;
            #endif

            // Verbose reporting of any detected error
            if (componentError) {
               std::complex<double> waveComponent =
                   std::complex<double>(dft_in[rank].x, dft_in[rank].y);
               Log::file() 
                  << " Error in FieldIo::convertKGridToBasis\n"
                  << "Star  Id        " << starPtr->starId << "\n"
                  << "Basis Id        " << starPtr->basisId << "\n"
                  << "Star WaveBz     " << starPtr->waveBz << "\n"
                  << "Wave indices    " << wavePtr->indicesDft << "\n"
                  << "Wave coeff      " << wavePtr->indicesDft << "\n"
                  << "Wave component  " << Dbl(waveComponent.real()) 
                                        << Dbl(waveComponent.imag()) << "\n"
                  << "Wave coefficient" << Dbl(coefficient.real()) 
                                        << Dbl(coefficient.imag()) << "\n"
                  << "Star component  " << Dbl(component.real()) 
                                        << Dbl(component.imag()) << "\n";
               UTIL_THROW("Error: Component error in FieldIo::convertKGridToBasis");
            }

            // Store real part
            components[ib] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either first wave of 1st star or last wave of 2nd star.
            wavePtr = &basis().wave(starPtr->beginId);
            if (wavePtr->implicit) {
               starPtr = &(basis().star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &basis().wave(starPtr->endId - 1);
               UTIL_CHECK(!(wavePtr->implicit));
               UTIL_CHECK(wavePtr->starId == is+1);
            }
            rank = dftMesh.rank(wavePtr->indicesDft);

            // Compute component value
            component =
                     std::complex<double>(dft_in[rank].x, dft_in[rank].y);
            UTIL_CHECK(abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);

            // Compute basis function coefficient values
            if (starPtr->invertFlag == 1) {
               components[ib] = component.real();
               components[ib+1] = -component.imag();
            } else {
               components[ib] = component.real();
               components[ib+1] = component.imag();
             }

            // Increment star counter by 2 (two stars were processed)
            is += 2;

         } else {

            UTIL_THROW("Invalid invertFlag value");

         }

      } //  loop over star index is

      // Deallocate arrays (clean up)
      delete[] dft_in;

   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray< DArray<double> > const & in,
                                        DArray< RDFieldDft<D> >& out) 
   const
   {
      UTIL_CHECK(in.capacity() == out.capacity());
      const int n = in.capacity();

      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }

   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(DArray< RDFieldDft<D> > const & in,
                                        DArray< DArray<double> > & out) 
   const
   {

      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();

      for (int i = 0; i < n; ++i) {
         convertKGridToBasis(in[i], out[i]);
      }

   }

   template <int D>
   void
   FieldIo<D>::convertBasisToRGrid(DArray< DArray<double> > const & in,
                                   DArray< RDField<D> >& out) const
   {
      UTIL_CHECK(basis().isInitialized());
      UTIL_CHECK(in.capacity() == out.capacity());
      const int nMonomer = in.capacity();

      DArray< RDFieldDft<D> > workDft;
      workDft.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         workDft[i].allocate(mesh().dimensions());
      }

      convertBasisToKGrid(in, workDft);

      for (int i = 0; i < nMonomer; ++i) {
         fft().inverseTransformSafe(workDft[i], out[i]);
         workDft[i].deallocate();
      }
   }

   template <int D>
   void FieldIo<D>::convertRGridToBasis(DArray< RDField<D> > const & in,
                                        DArray< DArray<double> > & out) 
   const
   {
      UTIL_CHECK(basis().isInitialized());
      UTIL_CHECK(in.capacity() == out.capacity());
      const int nMonomer = in.capacity();

      // Allocate work space
      DArray< RDFieldDft<D> > workDft;
      workDft.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         workDft[i].allocate(mesh().dimensions());
      }

      // Forward FFT, then KGrid to Basis
      for (int i = 0; i < nMonomer; ++i) {
         fft().forwardTransformSafe(in[i], workDft [i]);
      }
      convertKGridToBasis(workDft, out);

      // Deallocate work space
      for (int i = 0; i < nMonomer; ++i) {
         workDft[i].deallocate();
      }
   }

   template <int D>
   void
   FieldIo<D>::convertKGridToRGrid(DArray< RDFieldDft<D> > const & in,
                                   DArray< RDField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   template <int D>
   void
   FieldIo<D>::convertRGridToKGrid(DArray< RDField<D> > const & in,
                                   DArray< RDFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransformSafe(in[i], out[i]);
      }
   }

   template <int D>
   void FieldIo<D>::checkWorkDft() const
   {
      if (!workDft_.isAllocated()) {
         workDft_.allocate(mesh().dimensions());
      } else {
         UTIL_CHECK(workDft_.meshDimensions() == fft().meshDimensions());
      }
   }

} // namespace Pspc
} // namespace Pscf
#endif
