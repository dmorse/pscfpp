#ifndef RP_TRAJECTORY_WRITER_H
#define RP_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h>   // base class template

#include <iostream>
#include <fstream>

namespace Pscf {
namespace Rp {

   /**
   * Periodically write field frames (snapshots) to a trajectory file.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Cpp<D> or Rpg::Types<D>
   *
   * \see \ref rp_TrajectoryWriter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class TrajectoryWriter : public Analyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      TrajectoryWriter(Simulator<D,T>& simulator, 
                       System<D,T>& system);

      /**
      * Destructor.
      */
      ~TrajectoryWriter() = default;

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      void readParameters(std::istream& in) override;

      /**
      * Setup before main simulation loop. 
      */
      void setup() override;

      /**
      * Write a frame/snapshot to the trajectory file.
      *
      * \param iStep  step index
      */
      void sample(long iStep) override;

      /**
      * Close trajectory file after run.
      */
      void output() override;

   protected:

      /**
      * Write data that should appear once, at beginning of the file.
      *
      * \param out  output file stream
      */
      void writeHeader(std::ofstream& out);

      /**
      * Write data that should appear in every frame.
      *
      * \param out output file stream
      * \param iStep MC time step index
      */
      void writeFrame(std::ofstream& out, long iStep);

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

      // Inherited protected member functions (selected).
      using Analyzer<D,T>::simulator;
      using Analyzer<D,T>::system;
   };

}
}
#endif
