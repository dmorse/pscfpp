#ifndef RPC_CONCENTRATION_WRITER_H
#define RPC_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Periodically write snapshots to a trajectory file.
   *
   * \see \ref rpc_ConcentrationWriter_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationWriter : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      ConcentrationWriter(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~ConcentrationWriter()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear nSample counter.
      */
      virtual void setup();

      /**
      * Write a frame/snapshot to trajectory file.
      *
      * \param iStep step index
      */
      virtual void sample(long iStep);

      /**
      * Close trajectory file after run.
      */
      virtual void output();

      using ParamComposite::read;
      using ParamComposite::setClassName;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::isAtInterval;

   protected:

      // Output file stream
      std::ofstream outputFile_;

      // Output filename
      std::string filename_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;

      /// Has readParam been called?
      long isInitialized_;

      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

   protected:

      /**
      * Write data that should appear once, at beginning of the file.
      *
      * Called by sample on first invocation. Default implementation is empty.
      *
      * \param out output file stream
      */
      void writeHeader(std::ofstream& out);

      /**
      * Write data that should appear in every frame.
      *
      * \param out output file stream
      * \param iStep MC time step index
      */
      void writeFrame(std::ofstream& out, long iStep);

      /**
      * Return reference to parent system.
      */
      System<D>& system();

      /**
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();


   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& ConcentrationWriter<D>::system()
   {  return *systemPtr_; }

   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& ConcentrationWriter<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_CONCENTRATION_WRITER_TPP
   // Suppress implicit instantiation
   extern template class ConcentrationWriter<1>;
   extern template class ConcentrationWriter<2>;
   extern template class ConcentrationWriter<3>;
   #endif

}
}
#endif
