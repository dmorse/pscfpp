#ifndef PSPG_TRAJECTORY_WRITER_TPP
#define PSPG_TRAJECTORY_WRITER_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include "Analyzer.h"
#include <pspg/System.h>
#include <pspg/simulate/Simulator.h>
#include <util/misc/FileMaster.h>
//#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <sstream>


namespace Pscf {
namespace Rpg 
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

   
   #if 0
   /*
   * Load state from an archive.
   */
   template <int D>
   void TrajectoryWriter<D>::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      ar & nSample_;
      isInitialized_ = true;
   }
   #endif

   /*
   * Save state to archive.
   */
   template <int D>
   void TrajectoryWriter<D>::save(Serializable::OArchive& ar)
   { ar & *this; }

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
      UnitCell<D> const & unitCell = domain.unitCell();
      fieldIo.writeFieldsRGrid(out, system().w().rgrid(), unitCell,
                               writeHeader, isSymmetric);
      out << "\n";
   }  
   
   
    template <int D>
    void TrajectoryWriter<D>::writeHeader(std::ofstream& out)
   {  
      bool hasSymmetry = false;
      int nMonomer = system().mixture().nMonomer();
      FieldIo<D> const & fieldIo = system().domain().fieldIo(); 
      fieldIo.writeFieldHeader(out, nMonomer, system().domain().unitCell(), 
                               hasSymmetry);
      out << "\n";
                                       
   } 
   
   
   /*
   * Dump configuration to file
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
   * Read interval and outputFileName. 
   */
   template <int D>
   void TrajectoryWriter<D>::output() 
   {  outputFile_.close(); }

}
}
#endif
