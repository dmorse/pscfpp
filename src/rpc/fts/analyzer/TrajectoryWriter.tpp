#ifndef RPC_TRAJECTORY_WRITER_TPP
#define RPC_TRAJECTORY_WRITER_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include "Analyzer.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <util/misc/FileMaster.h>
//#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   TrajectoryWriter<D>::TrajectoryWriter(Simulator<D>& simulator, 
                                         System<D>& system) 
    : Analyzer<D>(),
      nSample_(0),
      isInitialized_(false),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {  setClassName("TrajectoryWriter"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void TrajectoryWriter<D>::readParameters(std::istream& in) 
   {
      Analyzer<D>::readParameters(in);
      isInitialized_ = true;
   }
   
   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void TrajectoryWriter<D>::setup() 
   {  
      nSample_ = 0; 
      std::string filename;
      filename_  = outputFileName();
      system().fileMaster().openOutputFile(filename_ , outputFile_);
      writeHeader(outputFile_);
   }

    template <int D>
    void TrajectoryWriter<D>::writeFrame(std::ofstream& out, long iStep)
   {  
      out << "i = " << iStep << "\n";
      bool writeHeader = false;
      bool isSymmetric = false;
      Domain<D> const & domain = system().domain();
      FieldIo<D> const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldsRGrid(out, system().w().rgrid(), 
                               domain.unitCell(), 
                               writeHeader, isSymmetric);
      out << "\n";
   }  
   
   
    template <int D>
    void TrajectoryWriter<D>::writeHeader(std::ofstream& out)
   {  
      int nMonomer = system().mixture().nMonomer();
      bool isSymmetric = false;
      Domain<D> const & domain = system().domain();
      FieldIo<D> const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldHeader(out, nMonomer, 
                               domain.unitCell(), isSymmetric);
      out << "\n";
   } 
   
   
   /*
   * Periodically write a frame to file
   */
   template <int D>
   void TrajectoryWriter<D>::sample(long iStep) 
   {  
      if (isAtInterval(iStep))  {
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }
  
   /*
   * Close output file at end of simulation.
   */
   template <int D>
   void TrajectoryWriter<D>::output() 
   {  outputFile_.close(); }

}
}
#endif
