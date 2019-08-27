/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <pscf/math/IntVec.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>
//#include <sstream>
//#include <unistd.h>

namespace Pscf {
namespace Pssp
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : mixturePtr_(0),
      unitCellPtr_(0),
      meshPtr_(0),
      groupNamePtr_(0),
      basisPtr_(0),
      fftPtr_(0),
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
   void FieldIo<D>::associate(Mixture<D>& mixture,
                             UnitCell<D>& unitCell,
                             Mesh<D>& mesh,
                             std::string& groupName,
                             Basis<D>& basis,
                             FFT<D>& fft,
                             FileMaster& fileMaster)
   {
      mixturePtr_ = &mixture;
      unitCellPtr_ = &unitCell;
      meshPtr_ = &mesh;
      groupNamePtr_ = &groupName;
      basisPtr_ = &basis;
      fftPtr_ = &fft;
      fileMasterPtr_ = &fileMaster;
   }

   template <int D>
   void FieldIo<D>::readFields(std::istream& in, 
                               DArray< DArray<double> >& fields)
   {
      readFieldHeader(in);

      std::string label;
      in >> label;
      UTIL_CHECK(label == "n_Star");
      int nStar;
      in >> nStar;
      UTIL_CHECK(nStar > 0);
      UTIL_CHECK(nStar == basis().nStar());

      // Read fields
      int nMonomer = mixture().nMonomer();
      int i, j, idum;
      for (i = 0; i < nStar; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nMonomer; ++j) {
            in >> fields[j][i];
         }
      }
   }
  
 
   template <int D>
   void FieldIo<D>::readFields(std::string filename, 
                              DArray<DArray<double> >& fields)
   {
       std::ifstream inFile;
       fileMaster().openInputFile(filename, inFile);
       readFields(inFile, fields);
       inFile.close();
   }


   template <int D>
   void FieldIo<D>::writeFields(std::ostream &out, 
                           DArray<DArray<double> > const &  fields)
   {
      writeFieldHeader(out);
      int nStar = basis().nStar();
      out << "nStar   "  <<  std::endl << "       " << nStar    <<std::endl;

      // Write fields
      int nMonomer = mixture().nMonomer();
      int i, j;
      for (i = 0; i < nStar; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nMonomer; ++j) {
            out << "  " << Dbl(fields[j][i], 18, 11);
         }
         //out<< "  " << basis().wave(basis().star(i).beginId).indicesDft;
         out << std::endl;
      }
   }


   template <int D>
   void FieldIo<D>::writeFields(std::string filename, 
                                DArray< DArray<double> > const &  fields)
   {
      std::ofstream outFile;
      fileMaster().openOutputFile(filename, outFile);
      writeFields(outFile, fields);
      outFile.close();
   }


   template <int D>
   void FieldIo<D>::readRFields(std::istream &in,
                                DArray<RField<D> >& fields)
   {
      readFieldHeader(in);

      std::string label;
      in >> label;
      UTIL_CHECK(label == "n_Grid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      int nMonomer = mixture().nMonomer();
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
      int n1 =0;
      int n2 =0;
      int n3 =0;

      if (D==3){
         while (n1 < mesh().dimension(0)){
            q = p;
            n2 = 0;
            while (n2 < mesh().dimension(1)){
               r =q;
               n3 = 0;
               while (n3 < mesh().dimension(2)){
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

      else if (D==2){
         while (n1 < mesh().dimension(0)){
            r =q; 
            n2 = 0;
            while (n2 < mesh().dimension(1)){
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

      else if (D==1){

         while (n1 < mesh().dimension(0)){
            for (int i = 0; i < nMonomer; ++i) {
               fields[i][s] = temp[i][r];
            }   
            ++r;
            ++s;
            ++n1;    
         }   

      } else {
         std::cout<<"Invalid Dimensions";
      }

   }


   template <int D>
   void FieldIo<D>::readKFields(std::istream &in,
                                DArray<RFieldDft<D> >& fields)
   {
      readFieldHeader(in);

      std::string label;
      in >> label;
      UTIL_CHECK(label == "n_Grid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());


      // Read Fields;
      int nMonomer = mixture().nMonomer();
      int idum;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         in >> idum;
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < 2; ++j) {
               in >> fields[i][itr.rank()][j];
            }
         }
      }
   }


   template <int D>
   void FieldIo<D>::writeRFields(std::ostream &out,
                                DArray<RField<D> > const& fields)
   {
      writeFieldHeader(out);
      out << "n_Grid" <<  std::endl
          << "            " << mesh().dimensions() << std::endl;

      int nMonomer = mixture().nMonomer();
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

      if (D==3){
         while (n3 < mesh().dimension(2)){
            q = p; 
            n2 = 0; 
            while (n2 < mesh().dimension(1)){
               r =q;
               n1 = 0; 
               while (n1 < mesh().dimension(0)){
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
      else if (D==2){
         while (n2 < mesh().dimension(1)){
            r =q;
            n1 = 0;
            while (n1 < mesh().dimension(0)){
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

      else if (D==1){

         while (n1 < mesh().dimension(0)){
            for (int i = 0; i < nMonomer; ++i) {
               temp[i][s] = fields[i][r];
            }
            ++r;
            ++s;
            ++n1;
         }
      }

      else{
         std::cout<<"Invalid Dimensions";
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
   void FieldIo<D>::writeKFields(std::ostream &out,
                           DArray<RFieldDft<D> > const& fields)
   {
      writeFieldHeader(out);
      out << "n_Grid" << std::endl 
          << "               " << mesh().dimensions() << std::endl;

      // Write fields
      int nMonomer = mixture().nMonomer();
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         out << Int(itr.rank(), 5);
         for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(fields[j][itr.rank()][0], 18, 11)
                   << Dbl(fields[j][itr.rank()][1], 18, 11);
         }
         out << std::endl;
      }
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
      UTIL_CHECK(nMonomer == mixture().nMonomer());
   }


   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out) const
   {
      out << "format  1   0" <<  std::endl;
      out << "dim" <<  std::endl 
          << "          " << D << std::endl;
      writeUnitCellHeader(out, unitCell()); 
      out << "group_name" << std::endl 
          << "          " << groupName() <<  std::endl;
      out << "N_monomer"  << std::endl 
          << "          " << mixture().nMonomer() << std::endl;
   }

   template <int D>
   void 
   FieldIo<D>::convertFieldBasisToDft(DArray<double> const& components, 
                                      RFieldDft<D>& dft)
   {
      // Create Mesh<D> with dimensions of DFT grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient of star
      std::complex<double> coeff;             // coefficient of wave
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int nStar = basis().nStar();            // number of stars
      int is;                                 // star index
      int iw;                                 // wave index

      is = 0;
      while (is < nStar) {
         starPtr = &(basis().star(is));
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Make real component (coefficient for star basis function)
            component = std::complex<double>(components[is], 0.0);

            // Loop over waves in closed star
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis().wave(iw));
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
            starPtr = &(basis().star(is));
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis().wave(iw));
               if (!(wavePtr->implicit)) {
                  coeff = component*(wavePtr->coeff);
                  indices = wavePtr->indicesDft;    
                  rank = dftMesh.rank(indices);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
               }
            }

            // Loop over waves in second star
            starPtr = &(basis().star(is+1));
            UTIL_CHECK(starPtr->invertFlag == -1);
            component = conj(component);
            for (iw = starPtr->beginId; iw < starPtr->endId; ++iw) {
               wavePtr = &(basis().wave(iw));
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
   void 
   FieldIo<D>::convertFieldDftToBasis(RFieldDft<D> const & dft, 
                                      DArray<double>& components)
   {
      // Create Mesh<D> with dimensions of DFT grid.
      Mesh<D> dftMesh(dft.dftDimensions());

      typename Basis<D>::Star const* starPtr; // pointer to current star
      typename Basis<D>::Wave const* wavePtr; // pointer to current wave
      std::complex<double> component;         // coefficient of star
      IntVec<D> indices;                      // dft grid indices of wave
      int rank;                               // dft grid rank of wave
      int nStar = basis().nStar();            // number of stars
      int is;                                 // star index

      // Loop over stars
      is = 0;
      while (is < nStar) {
         starPtr = &(basis().star(is));
         if (starPtr->cancel) continue;

         if (starPtr->invertFlag == 0) {

            // Characteristic wave is first wave of star
            wavePtr = &(basis().wave(starPtr->beginId));
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

} // namespace Pssp
} // namespace Pscf
