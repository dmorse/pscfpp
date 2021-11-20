#ifndef PSPC_FIELD_IO_TPP
#define PSPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : meshPtr_(0),
      fftPtr_(0),
      groupNamePtr_(0),
      basisPtr_(0),
      fileMasterPtr_()
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
                              std::string const & groupName,
                              Basis<D> const & basis,
                              FileMaster const & fileMaster)
   {
      meshPtr_ = &mesh;
      groupNamePtr_ = &groupName;
      basisPtr_ = &basis;
      fftPtr_ = &fft;
      fileMasterPtr_ = &fileMaster;
   }
  
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in, 
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Read header
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // Read nStar
      std::string label;
      in >> label;
      UTIL_CHECK(label == "N_star");
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      // Initialize all field array elements to zero
      int i, j;
      int nStarCapacity = fields[0].capacity();
      for (j = 0; j < nMonomer; ++j) {
         UTIL_CHECK(nStarCapacity == fields[j].capacity());
         for (i = 0; i < nStarCapacity; ++i) {
            fields[j][i] = 0.0;
         }
      }

      // Reset nStarIn = min(nStarIn, nStarCapacity)
      if (nStarCapacity < nStarIn) {
         nStarIn = nStarCapacity;
      }

      DArray<double> temp;
      temp.allocate(nMonomer);

      // Loop over stars to read field components
      IntVec<D> waveIn, waveBz, waveDft, waveStar;
      int waveId, starId, nWaveVectors;
      bool waveExists;
      for (i = 0; i < nStarIn; ++i) {

         // Read components for different monomers
         for (j = 0; j < nMonomer; ++j) {
            in >> temp [j];
         }

         // Read characteristic wave and number of wavectors in star.
         in >> waveIn;
         in >> nWaveVectors;

         // Check if waveIn is in first Brillouin zone (FBZ) for the mesh.
         waveBz = shiftToMinimum(waveIn, mesh().dimensions(), unitCell);
         waveExists = (waveIn == waveBz);

         // If wave is in FBZ, find in basis and set field components
         if (waveExists) {
            waveDft = waveBz;
            mesh().shift(waveDft);
            waveId = basis().waveId(waveDft);
            starId = basis().wave(waveId).starId;
            if (!basis().star(starId).cancel) {
               waveStar = basis().star(starId).waveBz;
               if (waveStar != waveBz) {
                   std::cout 
                     <<  "Inconsistent wave of star on input\n"
                     <<  "waveIn from file = " << waveIn   << "\n"
                     <<  "starId of waveIn = " << starId   << "\n"
                     <<  "waveBz of star   = " << waveStar << "\n";
                     UTIL_THROW("Inconsistent wave ids on file input");
               }
               UTIL_CHECK(basis().star(starId).waveBz == waveBz);
               for (j = 0; j < nMonomer; ++j) {
                  fields[j][starId] = temp [j];
               }
            }
         }

      }

   }
   
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::string filename, 
                                    DArray<DArray<double> >& fields,
                                    UnitCell<D>& unitCell)
   {
       std::ifstream file;
       fileMaster().openInputFile(filename, file);
       readFieldsBasis(file, fields, unitCell);
       file.close();
   }

   template <int D>
   void 
   FieldIo<D>::writeFieldsBasis(std::ostream &out, 
                                DArray<DArray<double> > const &  fields,
                                UnitCell<D> const & unitCell)
   const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Write header
      writeFieldHeader(out, nMonomer, unitCell);
      int nStar = basis().nStar();
      int nBasis = basis().nBasis();
      out << "N_star       " << std::endl 
          << "             " << nBasis << std::endl;

      // Write fields
      for (int i = 0; i < nStar; ++i) {
         if (!basis().star(i).cancel) {
            for (int j = 0; j < nMonomer; ++j) {
               out << Dbl(fields[j][i], 20, 10);
            }
            out << "   ";
            for (int j = 0; j < D; ++j) {
               out << Int(basis().star(i).waveBz[j], 5);
            } 
            out << Int(basis().star(i).size, 5) << std::endl;
         }
      }

   }

   template <int D>
   void 
   FieldIo<D>::writeFieldsBasis(std::string filename, 
                                DArray<DArray<double> > const & fields,
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
                                    DArray<RField<D> >& fields,
                                    UnitCell<D>& unitCell)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // Read grid dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp[i].allocate(mesh().dimensions());
      }

      // Read Fields;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int i = 0; i < nMonomer; ++i) {
            in  >> std::setprecision(15) >> temp[i][itr.rank()];
         }
      }

      int p = 0;
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;

      if (D==3) {
         while (n1 < mesh().dimension(0)) {
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               r = q;
               n3 = 0;
               while (n3 < mesh().dimension(2)) {
                  for (int i = 0; i < nMonomer; ++i) {
                     fields[i][s] = temp[i][r];
                  }
                  r = r + (mesh().dimension(0) * mesh().dimension(1));
                  ++s;
                  ++n3;              
               } 
               q = q + mesh().dimension(0);
               ++n2;
            } 
            ++n1;
            ++p;
         }
      }

      else if (D==2) {
         while (n1 < mesh().dimension(0)) {
            r =q; 
            n2 = 0;
            while (n2 < mesh().dimension(1)) {
               for (int i = 0; i < nMonomer; ++i) {
                  fields[i][s] = temp[i][r];
               }   
               r = r + (mesh().dimension(0));
               ++s;
               ++n2;    
            }   
            ++q;
            ++n1;
         }   
      } 

      else if (D==1) {

         while (n1 < mesh().dimension(0)) {
            for (int i = 0; i < nMonomer; ++i) {
               fields[i][s] = temp[i][r];
            }
            ++r;
            ++s;
            ++n1;    
         }   
      } 

      else{
         std::cout << "Invalid Dimensions";
      }

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename, 
                              DArray< RField<D> >& fields,
                              UnitCell<D>& unitCell)
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::ostream &out,
                                     DArray<RField<D> > const & fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      writeFieldHeader(out, nMonomer, unitCell);
      out << "ngrid" <<  std::endl
          << "           " << mesh().dimensions() << std::endl;

      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp[i].allocate(mesh().dimensions());
      } 

      int p = 0; 
      int q = 0; 
      int r = 0; 
      int s = 0; 
      int n1 =0;
      int n2 =0;
      int n3 =0;

      if (D==3) {
         while (n3 < mesh().dimension(2)) {
            q = p; 
            n2 = 0; 
            while (n2 < mesh().dimension(1)) {
               r =q;
               n1 = 0; 
               while (n1 < mesh().dimension(0)) {
                  for (int i = 0; i < nMonomer; ++i) {
                     temp[i][s] = fields[i][r];
                  }    
                  r = r + (mesh().dimension(1) * mesh().dimension(2));
                  ++s; 
                  ++n1;     
               }    
               q = q + mesh().dimension(2);
               ++n2;
            }    
            ++n3;
            ++p;     
         }    
      }
      else if (D==2) {
         while (n2 < mesh().dimension(1)) {
            r =q;
            n1 = 0;
            while (n1 < mesh().dimension(0)) {
               for (int i = 0; i < nMonomer; ++i) {
                  temp[i][s] = fields[i][r];
               }
               r = r + (mesh().dimension(1));
               ++s;
               ++n1;
            }
            ++q;
            ++n2;
         }
      }
      else if (D==1) {
         while (n1 < mesh().dimension(0)) {
            for (int i = 0; i < nMonomer; ++i) {
               temp[i][s] = fields[i][r];
            }
            ++r;
            ++s;
            ++n1;
         }
      } else {
         std::cout << "Invalid Dimensions";
      }

      // Write fields
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         // out << Int(itr.rank(), 5);
         for (int j = 0; j < nMonomer; ++j) {
            out << "  " << Dbl(temp[j][itr.rank()], 18, 15);
         }
         out << std::endl;
      }

   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::string filename, 
                                     DArray< RField<D> > const & fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::istream &in,
                                    DArray<RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   {
      // Inspect fields array
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
      }

      // Read header
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // Read mesh dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      // Read fields;
      int i, idum;
      MeshIterator<D> itr(fields[0].dftDimensions());
      i = 0;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         in >> idum;
         UTIL_CHECK(i == idum);
         UTIL_CHECK(i == itr.rank());
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < 2; ++j) {
               in >> fields[i][itr.rank()][j];
            }
         }
         ++i;
      }
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::string filename, 
                                    DArray< RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::ostream &out,
                                     DArray<RFieldDft<D> > const & fields,
                                     UnitCell<D> const & unitCell)
   const
   {
      // Inspect fields array
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
      }

      // Write header
      writeFieldHeader(out, nMonomer, unitCell);
      out << "ngrid" << std::endl 
          << "               " << mesh().dimensions() << std::endl;

      // Write fields
      MeshIterator<D> itr(fields[0].dftDimensions());
      int i = 0;
      for (itr.begin(); !itr.atEnd(); ++itr) {
         UTIL_CHECK(i == itr.rank());
         out << Int(itr.rank(), 5);
         for (int j = 0; j < nMonomer; ++j) {
               out << "  " 
                   << Dbl(fields[j][itr.rank()][0], 19, 12)
                   << Dbl(fields[j][itr.rank()][1], 19, 12);
         }
         out << std::endl;
         ++i;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::string filename, 
                                    DArray< RFieldDft<D> > const & fields,
                                    UnitCell<D> const & unitCell)
   const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   /*
   * Read common part of field header.
   */
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in, 
                                    int nMonomer,
                                    UnitCell<D>& unitCell) 
   {
      int ver1, ver2;
      std::string groupNameIn;
      int nMonomerIn;
      Pscf::readFieldHeader(in, ver1, ver2, unitCell, 
                            groupNameIn, nMonomerIn);
      // Note: Function definition in pscf/crystal/UnitCell.tpp
      if (groupNameIn != groupName()) {
         std::cout << std::endl 
             << "Warning - "
             << "Mismatched group names in FieldIo::readFieldHeader: \n" 
             << "  FieldIo::groupName :" << groupName() << "\n"
             << "  Field file header  :" << groupNameIn << "\n";
      }
      UTIL_CHECK(nMonomerIn == nMonomer);
   }

   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out, int nMonomer, 
                                     UnitCell<D> const & unitCell) const
   {
      int ver1 = 1;
      int ver2 = 0;
      Pscf::writeFieldHeader(out, ver1, ver2, unitCell, 
                             groupName(), nMonomer);
      // Note: This is defined in pscf/crystal/UnitCell.tpp
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray<double> const & in, 
                                        RFieldDft<D>& out)
   {
      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(out.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient for star
      std::complex<double> coeff;             // coefficient for wave
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int iw;                                 // wave index

      // Initialize all dft coponents to zero
      for (rank = 0; rank < dftMesh.size(); ++rank) {
         out[rank][0] = 0.0;
         out[rank][1] = 0.0;
      }

      // Loop over stars, skipping cancelled stars
      is = 0;
      while (is < basis().nStar()) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            component = std::complex<double>(in[is], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
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
            component = std::complex<double>(in[is], -in[is+1]);
            component /= sqrt(2.0);
            starPtr = &(basis().star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  out[rank][0] = coeff.real();
                  out[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis().star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = std::complex<double>(in[is], +in[is+1]);
            component /= sqrt(2.0);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
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

   template <int D>
   void FieldIo<D>::convertKGridToBasis(RFieldDft<D> const & in, 
                                        DArray<double>& out)
   {
      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(in.dftDimensions());

      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // coefficient for star
      int rank;                                // dft grid rank of wave
      int is;                                  // star index

      // Initialize all components to zero
      for (is = 0; is < basis().nStar(); ++is) {
         out[is] = 0.0;
      }

      // Loop over stars
      is = 0;
      while (is < basis().nStar()) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {
            ++is;
            continue;
         }

         if (starPtr->invertFlag == 0) {

            // Choose a wave in the star that is not implicit
            int beginId = starPtr->beginId;
            int endId = starPtr->endId;
            int iw = 0;
            bool isImplicit = true;
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
            component = std::complex<double>(in[rank][0], in[rank][1]);
            component /= wavePtr->coeff;
            UTIL_CHECK(abs(component.imag()) < 1.0E-8);
            out[is] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either the first wave of the 1st star or last wave of 2nd
            wavePtr = &basis().wave(starPtr->beginId);
            UTIL_CHECK(wavePtr->starId == is);
            if (wavePtr->implicit) {
               starPtr = &(basis().star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &basis().wave(starPtr->endId - 1);
               UTIL_CHECK(!(wavePtr->implicit));
               UTIL_CHECK(wavePtr->starId == is+1);
            }
            rank = dftMesh.rank(wavePtr->indicesDft);
            component = std::complex<double>(in[rank][0], in[rank][1]);
            UTIL_CHECK(abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);

            // Compute basis function coefficient values
            if (starPtr->invertFlag == 1) {
               out[is] = component.real();
               out[is+1] = -component.imag();
            } else {
               out[is] = component.real();
               out[is+1] = component.imag();
            }

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray< DArray <double> >& in,
                                        DArray< RFieldDft<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(DArray< RFieldDft<D> >& in,
                                        DArray< DArray <double> > & out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertKGridToBasis(in[i], out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertBasisToRGrid(DArray< DArray <double> >& in,
                                   DArray< RField<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkWorkDft();

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransform(workDft_, out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToBasis(DArray< RField<D> >& in,
                                   DArray< DArray <double> > & out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkWorkDft();

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], workDft_);
         convertKGridToBasis(workDft_, out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertKGridToRGrid(DArray< RFieldDft<D> >& in,
                                   DArray< RField<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransform(in[i], out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToKGrid(DArray< RField<D> >& in,
                                   DArray< RFieldDft<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Test if an RField<D> has declared space group symmetry.
   * Return true if symmetric, false otherwise.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RField<D> & in) 
   {
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_);
   }

   /*
   * Test if an RFieldDft has the declared space group symmetry.
   * Return true if symmetric, false otherwise.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RFieldDft<D> const & in) const
   {
      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> waveCoeff;         // coefficient from wave
      std::complex<double> rootCoeff;         // coefficient from root 
      std::complex<double> diff;              // coefficient difference
      int is;                                 // star index
      int iw;                                 // wave index
      int beginId, endId;                     // star begin, end ids
      int rank;                               // dft grid rank of wave

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(in.dftDimensions());

      // Loop over all stars
      for (is = 0; is < basis().nStar(); ++is) {
         starPtr = &(basis().star(is));

         if (starPtr->cancel) {

            // Check that coefficients are zero for all waves in star
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!wavePtr->implicit) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  if (abs(in[rank][0]) > 1.0E-9) return false;
                  if (abs(in[rank][1]) > 1.0E-9) return false;
               }
            }

         } else {

            // Check consistency of coeff values from all waves
            bool hasRoot = false;
            beginId = starPtr->beginId;
            endId = starPtr->endId;
            for (iw = beginId; iw < endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  rank = dftMesh.rank(wavePtr->indicesDft);
                  waveCoeff = std::complex<double>(in[rank][0], in[rank][1]);
                  waveCoeff /= wavePtr->coeff;
                  if (hasRoot) {
                     diff = waveCoeff - rootCoeff;
                     if (abs(diff) > 1.0E-9) return false;
                  } else {
                     rootCoeff = waveCoeff;
                     hasRoot = true;
                  }
               }
            }

         }

      } //  end loop over star index is

      // If the code reaches this point, the field is symmetric
      return true;
   }

   template <int D>
   void FieldIo<D>::checkWorkDft()
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
