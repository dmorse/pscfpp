#ifndef RPG_TRAJECTORY_READER_H
#define RPG_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    
#include <util/global.h>
#include <iostream>
#include <rpg/System.h>


namespace Pscf {
namespace Rpg
{

   template <int D> class System;
   using namespace Util;

   /**
   * Trajectory file reader.
   *
   * \ingroup Rpg_Fts_Trajectory_Module
   */
   template <int D>
   class TrajectoryReader
   {

   public:

      /**
      * Constructor.
      */
       TrajectoryReader<D>(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader(){};

      /**
      * Open trajectory file and read header, if any.
      *
      * By convention, this function treats the trajectory filename
      * as the name of an input file, and opens the file using the 
      * FileMaster:openInutFile function. This function prepends the 
      * input prefix (if any) to the file path. If compiled with MPI 
      * enabled, so that each processor simulates a different system, 
      * it also prepends a processor id prefix before the input prefix.
      *
      * \param filename trajectory input file name.
      */
      virtual void open(std::string filename) = 0;

      /**
      * Read a single frame. Frames are assumed to be read consecutively. 
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
      
      virtual void readHeader(){};
      
   protected:

      /** 
      * Return reference to parent system.
      */
      System<D>& system();
        
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;


   }; // end class TrajectoryReader
   

   // Get the parent system.
   template <int D>
   inline System<D>& TrajectoryReader<D>::system()
   {  return *systemPtr_; }
   
   
   

}
}
#endif
