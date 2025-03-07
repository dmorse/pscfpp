#ifndef RPG_ANALYZER_H
#define RPG_ANALYZER_H

#include <util/param/ParamComposite.h>      // base class
#include <util/misc/FileMaster.h>           // member variable

#include <string>
#include <iostream>
#include <fstream>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * The periodic action associated with an Analyzer may involve retrieval
   * or computation of a physical property, adding it to a statistical 
   * accumulator, and/or outputting it to file. This periodic action must
   * be implemented by the pure virtual sample() method.
   *
   * The sample() method should take the desired action only when the
   * simulation step index is an integer multiple of the associated interval
   * parameter.  The interval must be a positive integer that is a multiple 
   * of the static member Analyzer::baseInterval, which is set equal to 1 
   * by default.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class Analyzer : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Analyzer();

      /**
      * Destructor.
      */
      virtual ~Analyzer();

      /**
      * Read parameters from archive.
      *
      * Default implementation, reads interval and outputFileName.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      *
      * This method will be called just before the beginning of the
      * main simulation loop, after an initial field configuration 
      * is known. It may be used to complete any initialization that
      * cannot be completed in the readParameters method.
      *
      * The default implementation is an empty function.
      */
      virtual void setup()
      {}

      /**
      * Calculate, analyze and/or output a physical quantity.
      *
      * Take an action if iStep is a multiple of interval.
      * If iStep is not a multiple of interval, this method
      * should do nothing and return immediately.
      *
      * \param iStep current simulation step index.
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
      * Get interval value.
      */
      int interval() const;

      /**
      * Return true iff counter is a multiple of the interval.
      *
      * \param counter simulation step counter
      */
      bool isAtInterval(long counter) const;

      // Static members

      /**
      * The interval for every Analyzer must be a multiple of baseInterval.
      */
      static long baseInterval;

      /**
      * Define and initialize baseInterval.
      */
      static void initStatic();

   protected:

      /**
      * Set the FileMaster to use to open files.
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Read interval from file, with error checking.
      *
      * \param in input parameter file stream
      */
      void readInterval(std::istream &in);

      /**
      * Read outputFileName from file.
      *
      * \param in input parameter file stream.
      */
      void readOutputFileName(std::istream &in);

      /**
      * Get the FileMaster by reference.
      *
      * This can be used to open multiple output files.
      */
      FileMaster& fileMaster();

      /**
      * Return outputFileName string.
      */
      const std::string& outputFileName() const;

      /**
      * Return outputFileName string with added suffix.
      */
      std::string outputFileName(const std::string& suffix) const;

   private:

      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Base name of output file(s).
      std::string outputFileName_;
      
      /// Pointer to fileMaster for opening output file(s).
      FileMaster* fileMasterPtr_;

   };

   // Inline methods

   /*
   * Return interval value.
   */
   template <int D>
   inline int Analyzer<D>::interval() const
   {  return interval_; }

   /*
   * Return true iff the counter parameter is a multiple of the interval.
   */
   template <int D>
   inline bool Analyzer<D>::isAtInterval(long counter) const
   {  return (counter%interval_ == 0); }

   /*
   * Get the outputFileName string.
   */
   template <int D>
   inline const std::string& Analyzer<D>::outputFileName() const
   {  return outputFileName_; }

   // Method template

   #if 0
   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void Analyzer<D>::serialize(Archive& ar, const unsigned int version)
   {
      ar & interval_;
      ar & outputFileName_;
   }
   #endif

   #ifndef RPG_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class Analyzer<1>;
   extern template class Analyzer<2>;
   extern template class Analyzer<3>;
   #endif

}
}
#endif
