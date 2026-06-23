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
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named TrajectoryWriter, that
   * are defined in the Rpc and Rpg namespaces for use in the pscf_rpc and 
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_TrajectoryWriter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class TrajectoryWriter : public T::Analyzer
   {

   public:

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

      // Alias for base class.
      using AnalyzerT = typename T::Analyzer;

      // Inherited protected member functions (selected).
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

   };

}
}
#endif
