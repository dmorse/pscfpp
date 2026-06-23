#ifndef RP_ANALYZER_H
#define RP_ANALYZER_H

#include <util/param/ParamComposite.h>      // base class
#include <string>
#include <iostream>

// Forward declaration
namespace Util {
   class FileMaster;
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * Specializations of this class template are used as base classes for
   * two closely analogous classes, also named Analyzer, that are defined 
   * in the Rpc and Rpg namespaces for use in pscf_rpc and pscf_rpg.
   *
   * The periodic action associated with an Analyzer may involve retrieval
   * or computation of a physical property value, adding it to a statistical
   * accumulator, and/or outputting it to file. This periodic action must
   * be implemented by the pure virtual sample() function.
   *
   * The sample() function should take the desired action only when the
   * simulation step index is an integer multiple of the associated interval
   * parameter. The interval of each Analyzer must be a positive integer
   * that is a multiple of the static member Analyzer::baseInterval, which
   * is set to 1 by default.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - SimT : simulation type, e.g., Rpc::Simulation<D>
   *    - SysT : system type, e.g., Rpc::System<D>
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class SimT, class SysT>
   class Analyzer : public ParamComposite
   {

   public:

      // Protected constructor and destructor (see below).

      /**
      * Read parameters from archive.
      *
      * Default implementation, reads interval and outputFileName.
      *
      * \param in  input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      *
      * This function must be called just before the beginning of the main
      * simulation loop, after an initial field configuration is known.
      * It may be used to complete any initialization that cannot be
      * completed in the readParameters function, because knowledge of the
      * configuration is required.
      *
      * The default implementation is an empty function.
      */
      virtual void setup()
      {}

      /**
      * Calculate, analyze and/or output a physical quantity.
      *
      * Take an action if iStep is a multiple of the analyzer interval.
      * If iStep is not a multiple of interval, this function should do
      * nothing and return immediately.
      *
      * \param iStep  current simulation step index
      */
      virtual void sample(long iStep) = 0;

      /**
      * Output any results at the end of the simulation.
      *
      * The default implementation is an empty function.
      */
      virtual void output()
      {}

      /**
      * Get the interval value.
      */
      int interval() const;

      /**
      * Return true iff counter is a multiple of the interval.
      *
      * \param counter  simulation step counter
      */
      bool isAtInterval(long counter) const;

      // Static members

      /**
      * The interval for every Analyzer must be a multiple of baseInterval.
      */
      static long baseInterval;

      /**
      * Define and initialize baseInterval (initialized to 1).
      */
      static void initStatic();

   protected:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      Analyzer(SimT& simulator, SysT& system);

      /**
      * Destructor.
      */
      ~Analyzer() = default;

      /**
      * Set the FileMaster to use to open files.
      *
      * \param fileMaster  associated FileMaster object
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Optionally read interval from file, with error checking.
      *
      * If no interval parameter is present, the interval is set to 1
      * by default.
      *
      * \param in  input parameter file stream
      */
      void readInterval(std::istream& in);

      /**
      * Read outputFileName from file.
      *
      * \param in  input parameter file stream
      */
      void readOutputFileName(std::istream& in);

      /**
      * Get the parent Simulator by reference.
      */
      SimT& simulator();

      /**
      * Get the parent System by reference.
      */
      SysT& system();

      /**
      * Get the FileMaster by reference.
      *
      * This can be used to open multiple output files.
      */
      FileMaster& fileMaster();

      /**
      * Get the outputFileName string.
      */
      std::string const & outputFileName() const;

      /**
      * Return the outputFileName string with an added suffix.
      *
      * \param suffix  suffix that is appended to base outputFileName
      */
      std::string outputFileName(std::string const & suffix) const;

   private:

      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Base name of output file(s).
      std::string outputFileName_;

      /// Pointer to the parent Simulator.
      SimT* simulatorPtr_;

      /// Pointer to the parent System.
      SysT* systemPtr_;

      /// Pointer to fileMaster for opening output file(s).
      FileMaster* fileMasterPtr_;

   };

   // Inline member functions

   /*
   * Get the interval value.
   */
   template <int D, class SimT, class SysT> inline
   int Analyzer<D,SimT,SysT>::interval() const
   {  return interval_; }

   /*
   * Return true iff counter is a multiple of the interval.
   */
   template <int D, class SimT, class SysT> inline
   bool Analyzer<D,SimT,SysT>::isAtInterval(long counter) const
   {  return (counter%interval_ == 0); }

   /*
   * Get the outputFileName string.
   */
   template <int D, class SimT, class SysT> inline
   std::string const & Analyzer<D,SimT,SysT>::outputFileName() const
   {  return outputFileName_; }

   /*
   * Get the parent Simulator by reference.
   */
   template <int D, class SimT, class SysT> inline
   SimT& Analyzer<D,SimT,SysT>::simulator()
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   /*
   * Get the parent System by reference.
   */
   template <int D, class SimT, class SysT> inline
   SysT& Analyzer<D,SimT,SysT>::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

} // namespace Rp
} // namespace Pscf
#endif
