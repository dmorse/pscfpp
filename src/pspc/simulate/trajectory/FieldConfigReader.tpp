#ifndef PSPC_FIELD_CONFIG_READER_TPP
#define PSPC_FIELD_CONFIG_READER_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldConfigReader.h"
#include <pspc/System.h>
#include <pscf/math/IntVec.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/misc/ioUtil.h>
#include <iostream>
#include <string>

namespace Pscf {
namespace Pspc 
{

   using namespace Util;

   /*
   * Constructor. 
   */
   template <int D>
   FieldConfigReader<D>::FieldConfigReader(System<D>& system)
    : TrajectoryReader<D>(system),
      systemPtr_(&system),
      isAllocated_(false)
   {}

   template <int D>
   void FieldConfigReader<D>::allocate()
   {  
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      const int meshSize = system().domain().mesh().size();
      if (!isAllocated_){
         wField_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wField_[i].allocate(meshSize);
         }
         isAllocated_ = true;
      }
      meshDimensions_ = system().domain().mesh().dimensions();
   }
   
   /*
   * Open file and setup memory.
   */
   template <int D>
   void FieldConfigReader<D>::open(std::string filename)
   {
      system().fileMaster().open(filename, inputfile_);
      allocate();
   }
 
   template <int D>
   void FieldConfigReader<D>::readHeader()
   { 
      //Skip the header
      for (int i = 0; i < 13; ++i){
         inputfile_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      Log::file()<<"ReadHeader" << "\n";
      
   }
   /*
   * Read frame, return false if end-of-file
   */
   template <int D>
   bool FieldConfigReader<D>::readFrame()
   {
      // Preconditions
      if (!isAllocated_) {
         UTIL_THROW("Real Grid Field is not allocated");
      }
      
      bool notEnd;
      std::stringstream line;

      // Attempt to read first line
      notEnd = getNextLine(inputfile_, line);
      if (!notEnd) {
         return false;
      }
     
      // Process ITEM: TIMESTEP
      checkString(line, "i");
      checkString(line, "=");
      std::string value;
      line >> value;
      int step;
      step = std::stoi(value);
      Log::file()<< "step "<< step <<"\n";
      notEnd = getNextLine(inputfile_, line);
      checkString(line, "mesh");
      // Read ITEM: NUMBER OF Mesh
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: NUMBER OF Mesh");
      }
      notEnd = getNextLine(inputfile_, line);
      
      readFieldsRGrid(inputfile_, wField_);
      system().setWRGrid(wField_);

      return true;
   }
   
   /*
   * Close trajectory file.
   */
   template <int D>
   void FieldConfigReader<D>::close()
   {  inputfile_.close();}
   
   /*
   * read Frame of fields
   */
   
   template <int D>
   void FieldConfigReader<D>::readFieldsRGrid(std::istream &in,
                                              DArray<RField<D> >& fields)
   {
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> dimensions = system().domain().mesh().dimensions();
      // Setup temporary workspace array.
      DArray<RField<D> > temp;
      temp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         temp[i].allocate(dimensions);
      }

      // Read Fields;
      MeshIterator<D> itr(dimensions);
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
         while (n1 < system().domain().mesh().dimension(0)) {
            q = p;
            n2 = 0;
            while (n2 < system().domain().mesh().dimension(1)) {
               r = q;
               n3 = 0;
               while (n3 < system().domain().mesh().dimension(2)) {
                  for (int i = 0; i < nMonomer; ++i) {
                     fields[i][s] = temp[i][r];
                  }
                  r = r + (system().domain().mesh().dimension(0) * system().domain().mesh().dimension(1));
                  ++s;
                  ++n3;              
               } 
               q = q + system().domain().mesh().dimension(0);
               ++n2;
            } 
            ++n1;
            ++p;
         }
      }

      else if (D==2) {
         while (n1 < system().domain().mesh().dimension(0)) {
            r =q; 
            n2 = 0;
            while (n2 < system().domain().mesh().dimension(1)) {
               for (int i = 0; i < nMonomer; ++i) {
                  fields[i][s] = temp[i][r];
               }   
               r = r + (system().domain().mesh().dimension(0));
               ++s;
               ++n2;    
            }   
            ++q;
            ++n1;
         }   
      } 

      else if (D==1) {

         while (n1 < system().domain().mesh().dimension(0)) {
            for (int i = 0; i < nMonomer; ++i) {
               fields[i][s] = temp[i][r];
            }
            ++r;
            ++s;
            ++n1;    
         }   
      } 

   }


} 
}
#endif
