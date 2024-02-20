#ifndef PSPG_ANALYZER_TPP
#define PSPG_ANALYZER_TPP

#include "Analyzer.h"
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>


namespace Pscf {
namespace Rpg {

   using namespace Util;

   template <int D>
   void Analyzer<D>::initStatic()
   {  Analyzer<D>::baseInterval = 0; }

   /* 
   * Default constructor.
   */
   template <int D>
   Analyzer<D>::Analyzer()
    : ParamComposite(),
      outputFileName_(),
      interval_(1),
      fileMasterPtr_(0)
   {}
   
   /* 
   * Destructor.
   */
   template <int D>
   Analyzer<D>::~Analyzer()
   {}
   
   /*
   * Read parameters from stream, default implementation.
   */
   template <int D>
   void Analyzer<D>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
   }

   /*
   * Read the interval from parameter file, with error checking.
   */
   template <int D>
   void Analyzer<D>::readInterval(std::istream &in) 
   {
      // Check that baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }
   
      // Read interval value (inherited from Interval)
      read<long>(in, "interval", interval_);
   
      // Postconditons
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }
   }

   template <int D>
   void Analyzer<D>::readOutputFileName(std::istream &in) 
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   #if 0
   /*
   * Load parameters from archive, default implementation.
   */
   template <int D>
   void Analyzer<D>::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
   }

   /*
   * Load parameters from archive, with error checking.
   */
   template <int D>
   void Analyzer<D>::loadInterval(Serializable::IArchive& ar)
   {
      // Check that Analyzer<D>::baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }
   
      // Load parameters
      loadParameter<long>(ar, "interval", interval_);
   
      // Postconditons
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }
   }

   /*
   * Load outputFileName from archive.
   */
   template <int D>
   void Analyzer<D>::loadOutputFileName(Serializable::IArchive& ar)
   { loadParameter<std::string>(ar, "outputFileName", outputFileName_); }

   /*
   * Save interval and outputFileName to archive.
   */
   template <int D>
   void Analyzer<D>::save(Serializable::OArchive& ar)
   {
      ar & interval_;
      ar & outputFileName_;
   }
   #endif

   /*
   * Set the FileMaster.
   */
   template <int D>
   void Analyzer<D>::setFileMaster(FileMaster& fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Get the FileMaster by reference.
   */
   template <int D>
   FileMaster& Analyzer<D>::fileMaster()
   {  
      assert(fileMasterPtr_);
      return (*fileMasterPtr_);
   }

   /*
   * Get the outputFileName string with an added suffix
   */
   template <int D>
   std::string 
   Analyzer<D>::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

}
}
#endif 
