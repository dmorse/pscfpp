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
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : unitCellPtr_(0),
      meshPtr_(0),
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
   void FieldIo<D>::associate(UnitCell<D>& unitCell,
                             Mesh<D>& mesh,
                             FFT<D>& fft,
                             std::string& groupName,
                             Basis<D>& basis,
                             FileMaster& fileMaster)
   {
      unitCellPtr_ = &unitCell;
      meshPtr_ = &mesh;
      groupNamePtr_ = &groupName;
      basisPtr_ = &basis;
      fftPtr_ = &fft;
      fileMasterPtr_ = &fileMaster;
   }
  
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in, 
                                    DArray< RDField<D> >& fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      DArray<cudaReal*> temp_out;
      temp_out.allocate(nMonomer);

      // Read header
      FieldIo<D>::readFieldHeader(in);
      std::string label;
      in >> label;
      UTIL_CHECK(label == "N_star");
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      int i, j;
      int nStar = basis().nStar();
      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaReal[nStar];
      }   

      for (j = 0; j < nMonomer; ++j) {
         UTIL_CHECK(fields[j].capacity() == nStar);
         for (i = 0; i < nStar; ++i) {
            temp_out[j][i] = 0.0;
         }
      }

      DArray<double> temp;
      temp.allocate(nMonomer);

      // Loop over stars to read field components
      IntVec<D> waveIn, waveBz, waveDft;
      int waveId, starId, nWaveVectors;
      for (i = 0; i < nStarIn; ++i) {

         // Read components for different monomers
         for (j = 0; j < nMonomer; ++j) {
            in >> temp [j];
         }

         // Read characteristic wave and number of wavectors in star.
         in >> waveIn;
         in >> nWaveVectors;

         // Check if waveIn is in first Brillouin zone (FBZ) for the mesh.
         waveBz = shiftToMinimum(waveIn, mesh().dimensions(), unitCell());
         bool waveExists = (waveIn == waveBz);

         // If wave is in FBZ, find in basis and set field components
         if (waveExists) {
            waveDft = waveBz;
            mesh().shift(waveDft);
            waveId = basis().waveId(waveDft);
            starId = basis().wave(waveId).starId;
            UTIL_CHECK(basis().star(starId).waveBz == waveBz);
            if (!basis().star(starId).cancel) {
               for (j = 0; j < nMonomer; ++j) {
                  temp_out[j][starId] = temp [j];
               }
            }
         }

      }

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(fields[i].cDField(), temp_out[i],
            nStar * sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }

   }
   
 
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::string filename, 
                                    DArray<RDField<D> >& fields)
   {
       std::ifstream file;
       fileMaster().openInputFile(filename, file);
       readFieldsBasis(file, fields);
       file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsBasis(std::ostream &out, 
                                     DArray<RDField<D> > const &  fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      DArray<cudaReal*> temp_out;
      temp_out.allocate(nMonomer);

      // Write header
      writeFieldHeader(out, nMonomer);
      int nStar = basis().nStar();
      int nBasis = basis().nBasis();

      for(int i = 0; i < nMonomer; ++i) {
         temp_out[i] = new cudaReal[nStar];
      }   

     for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(temp_out[i], fields[i].cDField(),
            nStar * sizeof(cudaReal), cudaMemcpyDeviceToHost);
     }

      out << "N_star       " << std::endl 
          << "             " << nBasis << std::endl;

      // Write fields
      for (int i = 0; i < nStar; ++i) {
         if (!basis().star(i).cancel) {
            for (int j = 0; j < nMonomer; ++j) {
               out << Dbl(temp_out[j][i], 20, 10);
            }
            out << "   ";
            for (int j = 0; j < D; ++j) {
               out << Int(basis().star(i).waveBz[j], 5);
            } 
            out << Int(basis().star(i).size, 5) << std::endl;
         }
      }

     for(int i = 0; i < nMonomer; i++) {
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }   

   }

   template <int D>
   void FieldIo<D>::writeFieldsBasis(std::string filename, 
                                     DArray<RDField<D> > const & fields)
   {
       std::ofstream file;
       fileMaster().openOutputFile(filename, file);
       writeFieldsBasis(file, fields);
       file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::istream &in,
                                    DArray<RDField<D> >& fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      FieldIo<D>::readFieldHeader(in);

      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
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
      
      for(int i = 0; i < nMonomer; i++) {
         cudaMemcpy(fields[i].cDField(), temp_out[i],
            mesh().size() * sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp_out[i];
         temp_out[i] = nullptr;
      }

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename, 
                                    DArray< RDField<D> >& fields)
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::ostream &out,
                                     DArray<RDField<D> > const& fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      writeFieldHeader(out, nMonomer);
      out << "ngrid" <<  std::endl
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
                                     DArray< RDField<D> > const & fields)
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsRGrid(file, fields);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::istream &in, RDField<D> &field)
   {
      FieldIo<D>::readFieldHeader(in);

      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
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
      delete temp_out;
      temp_out = nullptr;

   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::string filename, RDField<D> &field)
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldRGrid(file, field);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::ostream &out, RDField<D> const & field)
   {

      writeFieldHeader(out, 1);
      out << "ngrid" <<  std::endl
          << "           " << mesh().dimensions() << std::endl;

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
      
      delete temp_out;
      temp_out = nullptr;
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::string filename, RDField<D> const & field)
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::istream &in,
                                    DArray<RDFieldDft<D> >& fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      DArray<cudaComplex*> temp_out;
      temp_out.allocate(nMonomer);      

      // Read header
      readFieldHeader(in);
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

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
                                    DArray< RDFieldDft<D> >& fields)
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::ostream &out,
                                     DArray<RDFieldDft<D> > const& fields)
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Write header
      writeFieldHeader(out, nMonomer);
      out << "ngrid" << std::endl 
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
                                     DArray< RDFieldDft<D> > const& fields)
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in) 
   {
      std::string label;

      in >> label;
      UTIL_CHECK(label == "format");
      int ver1, ver2;
      in >> ver1 >> ver2;
 
      in >> label;
      UTIL_CHECK(label == "dim");
      int dim;
      in >> dim;
      UTIL_CHECK(dim == D);

      readUnitCellHeader(in, unitCell());

      in >> label;
      UTIL_CHECK(label == "group_name");
      std::string groupName;
      in >> groupName;

      in >> label;
      UTIL_CHECK(label == "N_monomer");
      int nMonomer;
      in >> nMonomer;
      UTIL_CHECK(nMonomer > 0);
   }

   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out, int nMonomer) const
   {
      out << "format  1   0" <<  std::endl;
      out << "dim" <<  std::endl 
          << "          " << D << std::endl;
      writeUnitCellHeader(out, unitCell()); 
      out << "group_name" << std::endl 
          << "          " << groupName() <<  std::endl;
      out << "N_monomer"  << std::endl 
          << "          " << nMonomer << std::endl;
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(RDField<D> const& components, 
                                        RDFieldDft<D>& dft)
   {
     cudaReal* components_in;
     components_in = new cudaReal[basis().nStar()];
     cudaMemcpy(components_in, components.cDField(),
            basis().nStar() * sizeof(cudaReal), cudaMemcpyDeviceToHost);

     int kSize = 1;
     for (int i = 0; i < D; i++) {
        if (i == D - 1) {
           kSize *= (mesh().dimension(i) / 2 + 1); 
        }   
        else {
           kSize *= mesh().dimension(i);
        }   
     }   
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
      int iw;                                 // wave index

      // Initialize all dft coponents to zero
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

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            component = std::complex<double>(components_in[is], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!wavePtr->implicit) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = coeff.real();
                  dft_out[rank].y = coeff.imag();
               }
            }
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Make complex component for first star
            component = std::complex<double>(components_in[is], 
                                             -components_in[is+1]);
            component /= sqrt(2.0);

            // Loop over waves in first star
            starPtr = &(basis().star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = coeff.real();
                  dft_out[rank].y = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis().star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = conj(component);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &basis().wave(iw);
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;
                  rank = dftMesh.rank(indices);
                  dft_out[rank].x = coeff.real();
                  dft_out[rank].y = coeff.imag();
               }
            }

            // Increment is by 2 (two stars were processed)
            is += 2;

         } else {
 
            UTIL_THROW("Invalid invertFlag value");
  
         }
      }
    
     cudaMemcpy(dft.cDField(), dft_out,
              kSize * sizeof(cudaComplex), cudaMemcpyHostToDevice);
   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(RDFieldDft<D> const& dft, 
                                        RDField<D>& components)
   {
     cudaReal* components_out;
     components_out = new cudaReal[basis().nStar()];

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
     cudaMemcpy(dft_in, dft.cDField(),
            kSize * sizeof(cudaComplex), cudaMemcpyDeviceToHost);

      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // coefficient for star
      IntVec<D> indices;                       // dft grid indices of wave
      int rank;                                // dft grid rank of wave
      int is;                                  // star index
      int iw;                                  // wave id, within star 
      bool isImplicit;

      // Initialize all components to zero
      for (is = 0; is < basis().nStar(); ++is) {
         components_out[is] = 0.0;
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

            // Choose a characteristic wave that is not implicit.
            // Start with the first, alternately searching from
            // the beginning and end of star.
            isImplicit = true;
            iw = 0;
            while (isImplicit) {
                UTIL_CHECK(iw <= (starPtr->size)/2);
                wavePtr = &basis().wave(starPtr->beginId + iw);
                if (wavePtr->implicit) {
                   wavePtr = &basis().wave(starPtr->endId - 1 - iw);
                }
                isImplicit = wavePtr->implicit;
                ++iw;
            }
            UTIL_CHECK(wavePtr->starId == is);
            indices = wavePtr->indicesDft;
            rank = dftMesh.rank(indices);

            // Compute component value
            component = std::complex<double>(dft_in[rank].x, dft_in[rank].y);
            component /= wavePtr->coeff;
            UTIL_CHECK(abs(component.imag()) < 1.0E-7); // Did not satisfy 1.0E-8 accuracy due to float point accuracy
            components_out[is] = component.real();
            ++is;

         } else
         if (starPtr->invertFlag == 1) {

            // Identify a characteristic wave that is not implicit:
            // Either first wave of 1st star or last wave of 2nd star.
            wavePtr = &basis().wave(starPtr->beginId);
            if (wavePtr->implicit) {
               starPtr = &(basis().star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &basis().wave(starPtr->endId-1);
               UTIL_CHECK(!(wavePtr->implicit));
            } 
            indices = wavePtr->indicesDft;
            rank = dftMesh.rank(indices);

            // Compute component value
            component = std::complex<double>(dft_in[rank].x, dft_in[rank].y);
            UTIL_CHECK(abs(wavePtr->coeff) > 1.0E-8);
            component /= wavePtr->coeff;
            component *= sqrt(2.0);
            components_out [is] = component.real();
            components_out [is+1] = -component.imag();

            is += 2;
         } else {
            UTIL_THROW("Invalid invertFlag value");
         }

      } //  loop over star index is

     cudaMemcpy(components.cDField(), components_out,
            basis().nStar() * sizeof(cudaReal), cudaMemcpyHostToDevice);
   }

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray< RDField <D> >& in,
                                        DArray< RDFieldDft<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();

      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }

   }

   #if 0
   template <int D>
   void FieldIo<D>::convertKGridToBasis(RFieldDft<D> const & dft, 
                                      DArray<double>& components)
   {
      // Create Mesh<D> with dimensions of DFT grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient of star
      IntVec<D> indices;                      // dft grid indices of wave
      int nStar = basis().nStar();            // number of stars
      int rank;                               // dft grid rank of wave
      int is;                                 // star index
      int iw;                                 // wave id, within star
      bool isImplicit;

      // Loop over stars
      is = 0;
      while (is < nStar) {
         starPtr = &(basis().star(is));
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Choose a characteristic wave that is not implicit.
            // Start with the first, alternately searching from
            // the beginning and end of star.
            isImplicit = true;
            iw = 0;
            while (isImplicit) {
                UTIL_CHECK(iw <= (starPtr->size)/2);
                wavePtr = &basis().wave(starPtr->beginId + iw);
                if (wavePtr->implicit) {
                   wavePtr = &basis().wave(starPtr->endId - 1 - iw);
                }
                isImplicit = wavePtr->implicit;
                ++iw;
            }
            UTIL_CHECK(wavePtr->starId == is);
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
            wavePtr = &(basis().wave(starPtr->beginId));
            if (wavePtr->implicit) {
               starPtr = &(basis().star(is+1));
               UTIL_CHECK(starPtr->invertFlag == -1);
               wavePtr = &(basis().wave(starPtr->endId-1));
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
   #endif

   template <int D>
   void FieldIo<D>::convertKGridToBasis(DArray< RDFieldDft<D> >& in,
                                        DArray< RDField <D> > & out)
   {

      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();

      for (int i = 0; i < n; ++i) {
         convertKGridToBasis(in[i], out[i]);
      }   

   }

   template <int D>
   void 
   FieldIo<D>::convertBasisToRGrid(DArray< RDField<D> >& in,
                                   DArray< RDField<D> >& out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      //checkWorkDft();

      int nMonomer = in.capacity();

      DArray< RDFieldDft<D> > workDft;
      workDft.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         workDft[i].allocate(mesh().dimensions());
      } 

      convertBasisToKGrid(in, workDft);

      for (int i = 0; i < nMonomer; ++i) {
         fft().inverseTransform(workDft[i], out[i]);
         workDft[i].deallocate();
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToBasis(DArray< RDField<D> >& in,
                                   DArray< RDField<D> > & out)
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      //checkWorkDft();

      int nMonomer = in.capacity();

      DArray< RDFieldDft<D> > workDft;
      workDft.allocate(nMonomer);
      for(int i = 0; i < nMonomer; ++i) {
         workDft[i].allocate(mesh().dimensions());
      }   

      for (int i = 0; i < nMonomer; ++i) {
         fft().forwardTransform(in[i], workDft [i]);
      }

      convertKGridToBasis(workDft, out);

      for (int i = 0; i < nMonomer; ++i) {
         workDft[i].deallocate();
      }
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
