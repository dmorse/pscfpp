#ifndef RPG_FIELD_CONFIG_READER_TPP
#define RPG_FIELD_CONFIG_READER_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldConfigReader.h"
#include <rpg/System.h>
#include <pscf/math/IntVec.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/misc/ioUtil.h>
#include <iostream>
#include <string>

namespace Pscf {
namespace Rpg
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
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      if (!isAllocated_){
         wField_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wField_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
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
      // Read Header
      int nMonomer = system().mixture().nMonomer();
      Domain<D> & domain = system().domain();
      FieldIo<D> & fieldIo = domain.fieldIo();
      UnitCell<D> & unitcell = domain.unitCell();
      bool isSymmetric;
      fieldIo.readFieldHeader(inputfile_, nMonomer, unitcell, isSymmetric);
      Log::file() << "Read Header" << "\n";
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
      #if 0
      std::string value;
      line >> value;
      int step;
      step = std::stoi(value);
      Log::file()<< "step "<< step <<"\n";
      #endif
      
      // Read ITEM: NUMBER OF Mesh
      // notEnd = getNextLine(inputfile_, line);
      int nMonomer = system().mixture().nMonomer();
      // Read a single real-space grid field frame from trajectory file
      system().domain().fieldIo().readFieldsRGridData(inputfile_, wField_, nMonomer);
      // Update system real-space grid field 
      system().setWRGrid(wField_);

      return true;
   }
   
   /*
   * Close trajectory file.
   */
   template <int D>
   void FieldConfigReader<D>::close()
   {  inputfile_.close();}
   
} 
}
#endif
