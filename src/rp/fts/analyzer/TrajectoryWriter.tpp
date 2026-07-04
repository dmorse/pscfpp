#ifndef RP_TRAJECTORY_WRITER_TPP
#define RP_TRAJECTORY_WRITER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   TrajectoryWriter<D,T>::TrajectoryWriter(Simulator<D,T>& simulator,
                                         System<D,T>& system)
    : Analyzer<D,T>(simulator, system),
      nSample_(0),
      isInitialized_(false)
   {  ParamComposite::setClassName("TrajectoryWriter"); }

   /*
   * Read interval and outputFileName.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::readParameters(std::istream& in)
   {
      Analyzer<D,T>::readParameters(in);
      isInitialized_ = true;
   }

   /*
   * Setup before the main loop.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::setup()
   {
      nSample_ = 0;
      std::string filename = Analyzer<D,T>::outputFileName();
      system().fileMaster().openOutputFile(filename, outputFile_);
      writeHeader(outputFile_);
   }

   /*
   * Periodically write a frame to file.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::sample(long iStep)
   {
      if (Analyzer<D,T>::isAtInterval(iStep))  {
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }

   /*
   * Close the output file at end of simulation.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::output()
   {  outputFile_.close(); }

   /*
   * Write the trajectory file header.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::writeHeader(std::ofstream& out)
   {
      int nMonomer = system().mixture().nMonomer();
      bool isSymmetric = false;
      Domain<D,T> const & domain = system().domain();
      FieldIo<D,T> const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldHeader(out, nMonomer, domain.unitCell(),
                               isSymmetric);
      out << "\n";
   }

   /*
   * Write a frame in r-grid file format.
   */
   template <int D, class T>
   void TrajectoryWriter<D,T>::writeFrame(std::ofstream& out, long iStep)
   {
      out << "i = " << iStep << "\n";
      bool writeHeader = false;
      bool isSymmetric = false;
      Domain<D,T> const & domain = system().domain();
      FieldIo<D,T> const & fieldIo = domain.fieldIo();
      fieldIo.writeFieldsRGrid(out, system().w().rgrid(),
                               domain.unitCell(),
                               writeHeader, isSymmetric);
      out << "\n";
   }

}
}
#endif
