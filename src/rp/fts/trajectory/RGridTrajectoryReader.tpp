#ifndef RP_RGRID_TRAJECTORY_READER_TPP
#define RP_RGRID_TRAJECTORY_READER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridTrajectoryReader.h"
#include <prdc/crystal/UnitCell.h>
#include <util/misc/ioUtil.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Rp {

   template <int D, class T> class FieldIo;

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   RGridTrajectoryReader<D,T>::RGridTrajectoryReader(
                                         System<D,T>& system)
    : TrajectoryReader<D,T>(system),
      isAllocated_(false)
   {}

   /*
   * Allocate required memory.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::allocate()
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
   * Open file and allocate memory if necessary.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::open(std::string filename)
   {
      system().fileMaster().open(filename, inputfile_);
      allocate();
   }

   /*
   * Read trajectory file header section.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::readHeader()
   {
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D,T> const & fieldIo = system().domain().fieldIo();
      UnitCell<D> tmpUnitCell;
      bool hasSymmetry;
      fieldIo.readFieldHeader(inputfile_, nMonomer, tmpUnitCell,
                              hasSymmetry);
      system().setUnitCell(tmpUnitCell);
      Log::file() << "Read Header" << "\n";
   }

   /*
   * Read frame, return false if end-of-file.
   */
   template <int D, class T>
   bool RGridTrajectoryReader<D,T>::readFrame()
   {
      // Preconditions
      if (!isAllocated_) {
         UTIL_THROW("Memory not allocated in RGridTrajectoryReader");
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
      Log::file() << "step " << step << "\n";
      #endif

      // Read mesh dimensions
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);
      checkString(line, "mesh");

      // Read empty line
      notEnd = getNextLine(inputfile_, line);
      UTIL_CHECK(notEnd);

      // Read a w-field configuration in r-grid format
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D,T> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsRGridData(inputfile_, wField_, nMonomer);

      // Update system r-grid w fields
      system().w().setRGrid(wField_);

      return true;
   }

   /*
   * Close trajectory file.
   */
   template <int D, class T>
   void RGridTrajectoryReader<D,T>::close()
   {  inputfile_.close(); }

}
}
#endif
