#ifndef RP_RGRID_TRAJECTORY_READER_H
#define RP_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"           // base class
#include <pscf/math/IntVec.h>           // member
#include <util/containers/DArray.h>     // member
#include <fstream>                      // member
#include <string>                       // function parameter

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Trajectory file reader.
   *
   * Specializations of this class template are used as base classes for two
   * closely analogous class templates, also named RGridTrajectorReader,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \ingroup Rp_Fts_Trajectory_Module
   */
   template <int D, class T>
   class RGridTrajectoryReader : public T::TrajectoryReader
   {

   public:

      // Protected constructor and destructor (see below).

      /**
      * Open trajectory file and read header, if any.
      *
      * By convention, this function treats the trajectory filename
      * as the name of an input file, and opens the file using the
      * FileMaster:openInutFile function. This function prepends the
      * input prefix (if any) to the file path.
      *
      * \param filename  trajectory input file name
      */
      void open(std::string filename) override;

      /**
      * Read header of trajectory file (if any).
      */
      void readHeader() override;

      /**
      * Read a single frame from the trajectory file.
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      bool readFrame() override;

      /**
      * Close the trajectory file.
      */
      void close() override;

   protected:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      RGridTrajectoryReader(System<D,T>& system);

      /**
      * Destructor.
      */
      ~RGridTrajectoryReader() = default;

      using TrajectoryReaderT = typename T::TrajectoryReader;
      using TrajectoryReaderT::system;

   private:

      // Field configuration
      DArray< typename T::RField > wField_;

      // Dimensions of computational mesh (# of points in each direction)
      IntVec<D> meshDimensions_;

      // Trajectory file.
      std::ifstream inputfile_;

      // Has wField_ been allocated?
      bool isAllocated_;

      /**
      * Allocate required memory.
      */
      void allocate();

   };

}
}
#endif
