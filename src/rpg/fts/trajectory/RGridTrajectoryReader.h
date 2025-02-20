#ifndef RPG_RGRID_TRAJECTORY_READER_H
#define RPG_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"
#include <util/global.h>
#include <iostream>
#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <util/containers/DArray.h>

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
   class RGridTrajectoryReader : public TrajectoryReader<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System<D> object
      */
      RGridTrajectoryReader<D>(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~RGridTrajectoryReader(){};

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
      void open(std::string filename);

      /**
      * Read a single frame. Frames are assumed to be read consecutively. 
      *
      * This function reads a frame from the trajectory file that was
      * opened by the open() function.
      *
      * \return true if a frame is avaiable, false if at end of file
      */
      bool readFrame();

      /**
      * Close the trajectory file.
      */
      void close();
      
      /**
      * Read header of trajectory file.
      */
      void readHeader();
      
   protected:
      /**
      * Allocate memory required by trajectory reader
      */
      void allocate();
      
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      // Pointer to the parent system.
      System<D>* systemPtr_; 
      
      
   private:

      //
      IntVec<D> meshDimensions_;
      
      // Trajectory file.
      std::ifstream inputfile_;
      
      // Read Grid Field configuration      
      DArray< RField<D> > wField_;
      
      // Has the variable been allocated?
      bool isAllocated_;
        

   }; 
   
   // Get the parent system.
   template <int D>
   inline System<D>& RGridTrajectoryReader<D>::system()
   {  return *systemPtr_; }

   #ifndef RPG_RGRID_TRAJECTORY_READER_TPP
   // Suppress implicit instantiation
   extern template class RGridTrajectoryReader<1>;
   extern template class RGridTrajectoryReader<2>;
   extern template class RGridTrajectoryReader<3>;
   #endif
   
}
}
#endif
