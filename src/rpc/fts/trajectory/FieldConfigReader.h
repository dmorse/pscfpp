#ifndef RPC_FIELD_CONFIG_READER_H
#define RPC_FIELD_CONFIG_READER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"              // base class

#include <prdc/cpu/RField.h>               // member
#include <pscf/math/IntVec.h>              // member
#include <util/containers/DArray.h>        // member

#include <string>
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu; 

   /**
   * Trajectory file reader.
   *
   * \ingroup Rpc_Simulate_Trajectory_Module
   */
   template <int D>
   class FieldConfigReader : public TrajectoryReader<D>
   {

   public:

      /**
      * Constructor.
      */
      FieldConfigReader<D>(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~FieldConfigReader(){};

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
      
      void readFieldsRGrid(std::istream &in, DArray<RField<D> >& fields);

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
      
      // Readl Grid Field configuration      
      DArray< RField<D> > wField_;
      
      // Has the variable been allocated?
      bool isAllocated_;
        

   }; 
   
   // Get the parent system.
   template <int D>
   inline System<D>& FieldConfigReader<D>::system()
   {  return *systemPtr_; }

   
}
}
#endif
