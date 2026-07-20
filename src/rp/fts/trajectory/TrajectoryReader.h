#ifndef RP_TRAJECTORY_READER_H
#define RP_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
   }
}

namespace Pscf {
namespace Rp {

   /**
   * Trajectory file reader (abstract base class).
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \ingroup Rp_Fts_Trajectory_Module
   */
   template <int D, class T>
   class TrajectoryReader
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      TrajectoryReader(System<D,T>& system)
       : systemPtr_(&system)
      {}

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader() = default;

      /**
      * Open trajectory file and allocate memory if necessary.
      *
      * By convention, this function treats the trajectory filename
      * as the name of an input file, and opens the file using the
      * FileMaster:openInutFile function. This function prepends the
      * input prefix (if any) to the file path.
      *
      * \param filename trajectory input file name.
      */
      virtual void open(std::string filename) = 0;

      /**
      * Read header of trajectory file (if any).
      *
      * Empty default implementation.
      */
      virtual void readHeader()
      {};

      /**
      * Read a single frame.
      *
      * This function reads a frame from the trajectory file that was
      * opened by the open() function.
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      virtual bool readFrame() = 0;

      /**
      * Close the trajectory file.
      */
      virtual void close() = 0;

   protected:

      /**
      * Return reference to parent system.
      */
      System<D,T>& system()
      {  return *systemPtr_; }

   private:

      /**
      * Pointer to the parent system (not owned by this).
      */
      System<D,T>* systemPtr_;

   };

}
}
#endif
