#ifndef RPC_RGRID_TRAJECTORY_READER_TPP
#define RPC_RGRID_TRAJECTORY_READER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridTrajectoryReader.h"

#include <rpc/System.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/misc/ioUtil.h>

#include <sstream>
#include <iostream>
#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RGridTrajectoryReader<D>::RGridTrajectoryReader(System<D>& system)
    : TrajectoryReader<D>(system),
      isAllocated_(false)
   {}

   template <int D>
   void RGridTrajectoryReader<D>::allocate()
   {
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);

      meshDimensions_ = system().domain().mesh().dimensions();
      if (!isAllocated_){
         wField_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wField_[i].allocate(meshDimensions_);
         }
         isAllocated_ = true;
      }
   }

   /*
   * Open file and setup memory.
   */
   template <int D>
   void RGridTrajectoryReader<D>::open(std::string filename)
   {
      system().fileMaster().open(filename, inputfile_);
      allocate();
   }

   template <int D>
   void RGridTrajectoryReader<D>::readHeader()
   {
      // Read Header
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      UnitCell<D> tmpUnitCell;
      bool hasSymmetry;
      fieldIo.readFieldHeader(inputfile_, nMonomer, tmpUnitCell,
                              hasSymmetry);
      system().setUnitCell(tmpUnitCell);
      Log::file() << "Read Header" << "\n";
   }

   /*
   * Read frame, return false if end-of-file
   */
   template <int D>
   bool RGridTrajectoryReader<D>::readFrame()
   {
      // Preconditions
      if (!isAllocated_) {
         UTIL_THROW("R-grid Field is not allocated");
      }

      bool notEnd;
      std::stringstream line;

      // Attempt to read first line, check for end of file
      notEnd = getNextLine(inputfile_, line);
      if (!notEnd) {
         return false;
      }

      // Read line containing time step
      checkString(line, "i");
      checkString(line, "=");
      #if 0
      std::string value;
      line >> value;
      int step;
      step = std::stoi(value);
      Log::file()<< "step "<< step <<"\n";
      #endif

      // Read mesh dimensions 
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);
      checkString(line, "mesh");

      // Read empty line
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);

      // Read w-field configuration in r-grid format
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGridData(inputfile_, wField_, nMonomer);

      // Update system r-grid field
      system().setWRGrid(wField_);

      return true;
   }

   /*
   * Close trajectory file.
   */
   template <int D>
   void RGridTrajectoryReader<D>::close()
   {  inputfile_.close();}

}
}
#endif
