#ifndef PSPC_FIELD_IO_TPP
#define PSPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>

#include <util/misc/Log.h>
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
   void FieldIo<D>::readFieldBasis(std::istream& in, DArray<double>& field,
                                   UnitCell<D>& unitCell) 
   const
   {
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());
      readFieldsBasis(in, fields, unitCell);
      field = fields[0];
   }

   template <int D>
   void FieldIo<D>::readFieldBasis(std::string filename, 
                                   DArray<double>& field,
                                   UnitCell<D>& unitCell) 
   const
   {
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());
      readFieldsBasis(filename, fields, unitCell);
      field = fields[0];
   }

   template <int D>
   void FieldIo<D>::writeFieldBasis(std::ostream& out, 
                                    DArray<double> const & field,
                                    UnitCell<D> const & unitCell) 
   const
   {
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());
      fields[0] = field;
      writeFieldsBasis(out, fields, unitCell);
   }

   template <int D>
   void FieldIo<D>::writeFieldBasis(std::string filename, 
                                    DArray<double> const & field,
                                    UnitCell<D> const & unitCell) 
   const
   {
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());
      fields[0] = field;
      writeFieldsBasis(filename, fields, unitCell);
   }

   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in, 
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell) 
   const
   {

      // Read in the field header and get the "number of monomers", 
      // equivalent to the number of fields in the file.
      int nMonomer;
      int fieldCapacity;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // Read the number of stars into nStarIn
      std::string label;
      in >> label;
      UTIL_CHECK(label == "N_star");
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      // If "fields" passed by reference is already allocated, check that the
      // input stream field dimensions match the fields parameter. 
      
      // Otherwise, allocate fields parameter to match the input stream.
      
      if (fields.isAllocated()) {
         // If outer DArray is allocated, require that it matches the 
         // number of inputted fields and that internal DArrays are also
         // allocated all with the same dimensions.
         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);

         // Check that all internal DArrays have same dimension.
         fieldCapacity = fields[0].capacity();
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK( fields[i].capacity() == fieldCapacity );
         }
      } else {
         // Else, allocate fields to the number of inputted fields
         // and the internal dimensions to the number of stars.
         fields.allocate(nMonomer);
         fieldCapacity = nStarIn;

         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(fieldCapacity);
         }
      }    

      // Initialize all field array elements to zero
      for (int i = 0; i < nMonomer; ++i) {
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
      int starId, starId2;
      int basisId, basisId2;
      int waveId, waveId2;

      std::complex<double> coeff, phasor;
      IntVec<D> waveBz, waveDft;
      int nWaveVector;
      int nReversedPair = 0;
      bool waveExists;

      // Loop over stars in input file to read field components
      int i = 0;
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
            //basisId = starId;
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

            } else {

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
                                    DArray<DArray<double> >& fields,
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
      int ib = 0;
      for (int i = 0; i < nStar; ++i) {
         if (!basis().star(i).cancel) {
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
   const
   {
      // Read in the field header and get the "number of monomers", 
      // equivalent to the number of fields in the file.
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);


      // If "fields" passed by reference is already allocated, check it is allocated
      // to match the system's mesh. 
      
      // Otherwise, allocate fields parameter to match the system's mesh.
      
      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);

         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(mesh().dimensions());
         }
      }

      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      // Setup temporary workspace array.
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
         Log::file() << "Invalid Dimensions";
      }

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename, 
                              DArray< RField<D> >& fields,
                              UnitCell<D>& unitCell)
   const
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
         Log::file() << "Invalid Dimensions";
      }

      // Write fields
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
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
   void FieldIo<D>::readFieldRGrid(std::istream &in, 
                                    RField<D> & field, 
                                    UnitCell<D>& unitCell)
   const
   {
      // Read in the field header and get the "number of monomers", 
      // equivalent to the number of fields in the file.
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // Only reading in a file with a single field.
      UTIL_CHECK(nMonomer == 1);

      // If "fields" passed by reference is already allocated, check it is allocated
      // to match the system's mesh. 
      
      // Otherwise, allocate fields parameter to match the system's mesh.

      if (field.isAllocated()) {
         UTIL_CHECK(field.meshDimensions() == mesh().dimensions());
      } else {
         field.allocate(mesh().dimensions());
      }

      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      // Setup temporary workspace.
      RField<D> temp;
      temp.allocate(mesh().dimensions());

      // Read Field;
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         in  >> std::setprecision(15) >> temp[itr.rank()];
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
                  field[s] = temp[r];
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
               field[s] = temp[r];
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
               field[s] = temp[r];
            ++r;
            ++s;
            ++n1;    
         }   
      } 

      else{
         Log::file() << "Invalid Dimensions";
      }
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::string filename, 
                                    RField<D> & field, 
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
                                    RField<D> const & field, 
                                    UnitCell<D> const & unitCell,
                                    bool writeHeader)
   const
   {
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell);
         out << "ngrid" <<  std::endl
             << "           " << mesh().dimensions() << std::endl;
      }

      RField<D> temp;
      temp.allocate(mesh().dimensions());

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
                  temp[s] = field[r];
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
               temp[s] = field[r];
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
            temp[s] = field[r];
            ++r;
            ++s;
            ++n1;
         }
      } else {
         Log::file() << "Invalid Dimensions";
      }

      // Write field
      MeshIterator<D> itr(mesh().dimensions());
      for (itr.begin(); !itr.atEnd(); ++itr) {
         out << "  " << Dbl(temp[itr.rank()], 18, 15);
         out << std::endl;
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::string filename, 
                                    RField<D> const & field, 
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
                                    DArray<RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      // Read in the field header and get the "number of monomers", 
      // equivalent to the number of fields in the file.
      int nMonomer;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell);

      // If "fields" passed by reference is already allocated, check it is allocated
      // to match the system's mesh. 
      
      // Otherwise, allocate fields parameter to match the system's mesh.
      
      if (fields.isAllocated()) {

         int nMonomerFields = fields.capacity();

         UTIL_CHECK(nMonomerFields > 0)
         UTIL_CHECK(nMonomerFields == nMonomer)

         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
         }
         
      } else {

         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(mesh().dimensions());
         }

      }      

      // Read and check input stream mesh dimensions
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
   const
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
                   << Dbl(fields[j][itr.rank()][0], 20, 12)
                   << Dbl(fields[j][itr.rank()][1], 20, 12);
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
   * Read common part of field header and extract 
   * the number of monomers (number of fields) in the file.
   */
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in, 
                                    int& nMonomer,
                                    UnitCell<D>& unitCell) 
   const
   {
      int ver1, ver2;
      std::string groupNameIn;
      UnitCell<D> tempUnitCell;
      tempUnitCell = unitCell; // make duplicate unit cell
      Pscf::readFieldHeader(in, ver1, ver2, unitCell, 
                            groupNameIn, nMonomer);
      // Note: Function definition in pscf/crystal/UnitCell.tpp

      // if the unit cell that was passed into this function was set 
      // (if nParameter > 0), then check that the field header data
      // matches the data that was originally in the unitCell
      if (tempUnitCell.nParameter() > 0) {

         // Check whether crystal system matches expectation
         if (unitCell.lattice() != tempUnitCell.lattice()) {
            Log::file() << std::endl 
               << "Mismatched crystal systems in FieldIo::readFieldHeader: \n" 
               << "  Expected crystal system : " << tempUnitCell.lattice() 
               << "\n  Field file header       : " << unitCell.lattice() 
               << std::endl;
            UTIL_THROW("Mismatched crystal system in field file header");
         }

         // Check whether lattice parameters match expectation
         if (unitCell.nParameter() != tempUnitCell.nParameter()) {
            // Check if nParameters is matched
            Log::file() << std::endl 
               << "Mismatched number of lattice parameters in "
               << "FieldIo::readFieldHeader: \n" 
               << "  Expected N_cell_param : " << tempUnitCell.nParameter()
               << "\n  Field file header     : " << unitCell.nParameter() 
               << std::endl;
            UTIL_THROW("Mismatched N_cell_param value in field file header");
         } else {
            // Print notice that lattice parameters are being overwritten
            Log::file() << std::endl
               << "Using lattice parameters from field file header.\n"
               << "Discarding previous lattice parameters." << std::endl;
         }

         // Check whether space group matches expectation
         if (groupNameIn != groupName()) {
            Log::file() << std::endl 
               << "Mismatched group names in FieldIo::readFieldHeader: \n" 
               << "  Expected group name : " << groupName() << "\n"
               << "  Field file header   : " << groupNameIn << "\n"
               << std::endl;
            UTIL_THROW("Mismatched space group name in field file header");
         }

      }
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
                                        RFieldDft<D>& out) const
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
      int ib;                                 // basis index
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

         // Set basisId for uncancelled star
         ib = starPtr->basisId;

         if (starPtr->invertFlag == 0) {

            // Make complex coefficient for star basis function
            component = std::complex<double>(in[ib], 0.0);

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
            component = std::complex<double>(in[ib], -in[ib+1]);
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
            component = std::complex<double>(in[ib], +in[ib+1]);
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
                                        DArray<double>& out) const
   {
      // Create Mesh<D> with dimensions of DFT Fourier grid.
      Mesh<D> dftMesh(in.dftDimensions());

      typename Basis<D>::Star const* starPtr;  // pointer to current star
      typename Basis<D>::Wave const* wavePtr;  // pointer to current wave
      std::complex<double> component;          // coefficient for star
      int rank;                                // dft grid rank of wave
      int is;                                  // star index
      int ib;                                  // basis index

      // Check if kgrid has symmetry, and print a warning if it does not
      bool symmetric;
      symmetric = hasSymmetry(in, true);
      if (!symmetric) {
         Log::file() << "WARNING: non-negligible error in conversion to "
                     << "symmetry-adapted basis format." 
                     << std::endl
                     << "See error values printed above."
                     << "The field output by this operation will be "
                     << std::endl
                     << "a symmetrized version of the input field."
                     << std::endl;
      }

      // Initialize all components to zero 
      for (is = 0; is < basis().nBasis(); ++is) {
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

         // Set basis id for uncancelled star
         ib = starPtr->basisId;

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
            out[ib] = component.real();
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

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray< DArray <double> > const & in,
                                        DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(DArray< RFieldDft<D> > const & in,
                                        DArray< DArray <double> > & out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertKGridToBasis(in[i], out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertBasisToRGrid(DArray<double> const & in, 
                                   RField<D>& out) const
   {
      checkWorkDft();
      convertBasisToKGrid(in, workDft_);
      fft().inverseTransformSafe(workDft_, out);
   }

   template <int D>
   void 
   FieldIo<D>::convertBasisToRGrid(DArray< DArray <double> > const & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkWorkDft();

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransformSafe(workDft_, out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToBasis(RField<D> const & in,
                                   DArray<double> & out) const
   {
      checkWorkDft();
      fft().forwardTransform(in, workDft_);
      convertKGridToBasis(workDft_, out);
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToBasis(DArray< RField<D> > const & in,
                                   DArray< DArray <double> > & out) const
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
   FieldIo<D>::convertKGridToRGrid(DArray< RFieldDft<D> > & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertKGridToRGrid(RFieldDft<D>& in, RField<D>& out) const
   {
      fft().inverseTransformSafe(in, out);
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToKGrid(DArray< RField<D> > const & in,
                                   DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   template <int D>
   void 
   FieldIo<D>::convertRGridToKGrid(RField<D> const & in,
                                   RFieldDft<D>& out) const
   {
      fft().forwardTransform(in, out);
   }

   /*
   * Test if an RField<D> has declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RField<D> const & in, bool verbose) const
   {
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_, verbose);
   }

   /*
   * Test if an RFieldDft has the declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RFieldDft<D> const & in, bool verbose) const
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

      double cancelledError(0.0);   // max error from cancelled stars
      double uncancelledError(0.0); // max error from uncancelled stars

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
                  waveCoeff = std::complex<double>(in[rank][0], in[rank][1]);
                  if (std::abs(waveCoeff) > cancelledError) {
                     cancelledError = std::abs(waveCoeff);
                     if ((!verbose) && (cancelledError > 1e-8)) return false;
                  }
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
                     if (std::abs(diff) > uncancelledError) {
                        uncancelledError = std::abs(diff);
                        if ((!verbose) && (uncancelledError > 1e-8)) {
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

      if ((cancelledError < 1e-8) && (uncancelledError < 1e-8)) {
         return true;
      } else if (verbose) {
         Log::file() << "Maximum coefficient of a cancelled star: "
                     << cancelledError << std::endl
                     << "Maximum error of coefficient for uncancelled star: "
                     << uncancelledError << std::endl;
      }
      return false;
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
