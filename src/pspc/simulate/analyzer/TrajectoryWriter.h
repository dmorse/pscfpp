#ifndef PSPC_TRAJECTORY_WRITER_H
#define PSPC_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Periodically write snapshots to a trajectory file.
   *
   * \ingroup Pspc_Simulate_Analyzer_Module
   */
   template <int D>
   class TrajectoryWriter : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      TrajectoryWriter(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~TrajectoryWriter()
      {}

      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      #if 0
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive.
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
      #endif

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

   #if 0
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void TrajectoryWriter<D>::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSample_;
   }
   #endif

   // Get the parent system.
   template <int D>
   inline System<D>& TrajectoryWriter<D>::system()
   {  return *systemPtr_; }

   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& TrajectoryWriter<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef PSPC_TRAJECTORY_WRITER_TPP
   // Suppress implicit instantiation
   extern template class TrajectoryWriter<1>;
   extern template class TrajectoryWriter<2>;
   extern template class TrajectoryWriter<3>;
   #endif

}
}
#endif
