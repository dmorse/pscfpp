#ifndef RP_CONCENTRATION_WRITER_H
#define RP_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h> // base class template

#include <pscf/backend/TmplDeclare.h>
#include <iostream>
#include <fstream>

namespace Pscf {
namespace Rp {

   /**
   * Periodically write c-field snapshots to a trajectory file.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_ConcentrationWriter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class ConcentrationWriter : public Analyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationWriter(Simulator<D,T>& simulator, 
                          System<D,T>& system);

      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize before main simulation loop.
      */
      void setup() override;

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep  step index
      */
      void sample(long iStep) override;

      /**
      * Close trajectory file after run.
      */
      void output() override;

      using AnalyzerT = Analyzer<D,T>;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   private:

      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

      /**
      * Write data that should appear once, at beginning of the file.
      *
      * \param out  output file stream
      */
      void writeHeader(std::ofstream& out);

      /**
      * Write data that should appear in every frame.
      *
      * \param out  output file stream
      * \param iStep  step index
      */
      void writeFrame(std::ofstream& out, long iStep);

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(ConcentrationWriter)

}
}
#endif
