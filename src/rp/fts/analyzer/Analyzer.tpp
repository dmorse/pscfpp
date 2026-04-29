#ifndef RP_ANALYZER_TPP
#define RP_ANALYZER_TPP

#include "Analyzer.h"
#include <util/misc/FileMaster.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Static members

   template <int D, class SimT, class SysT>
   long Analyzer<D,SimT,SysT>::baseInterval = 1;

   /*
   * Static initialization function.
   */
   template <int D, class SimT, class SysT>
   void Analyzer<D,SimT,SysT>::initStatic()
   {  Analyzer<D,SimT,SysT>::baseInterval = 1; }

   // Non-static member functions

   /*
   * Constructor.
   */
   template <int D, class SimT, class SysT>
   Analyzer<D,SimT,SysT>::Analyzer(SimT& simulator, SysT& system)
    : ParamComposite(),
      interval_(1),
      outputFileName_(""),
      simulatorPtr_(&simulator),
      systemPtr_(&system),
      fileMasterPtr_(nullptr)
   {}

   /*
   * Read parameters from stream, default implementation.
   */
   template <int D, class SimT, class SysT>
   void Analyzer<D,SimT,SysT>::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
   }

   /*
   * Read the interval from parameter file, with error checking.
   */
   template <int D, class SimT, class SysT>
   void Analyzer<D,SimT,SysT>::readInterval(std::istream& in)
   {
      // Check that baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }

      // Optionally read interval value (set to 1 by default)
      interval_ = 1;
      ParamComposite::readOptional<long>(in, "interval", interval_);

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

   template <int D, class SimT, class SysT>
   void Analyzer<D,SimT,SysT>::readOutputFileName(std::istream &in)
   {
      ParamComposite::read<std::string>(in, "outputFileName",
                                        outputFileName_);
   }

   /*
   * Set the FileMaster.
   */
   template <int D, class SimT, class SysT>
   void Analyzer<D,SimT,SysT>::setFileMaster(FileMaster& fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Get the FileMaster by reference.
   */
   template <int D, class SimT, class SysT>
   FileMaster& Analyzer<D,SimT,SysT>::fileMaster()
   {
      UTIL_CHECK(fileMasterPtr_);
      return (*fileMasterPtr_);
   }

   /*
   * Get the outputFileName string with an added suffix
   */
   template <int D, class SimT, class SysT>
   std::string
   Analyzer<D,SimT,SysT>::outputFileName(std::string const & suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

}
}
#endif
