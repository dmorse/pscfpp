#ifndef RPC_ANALYZER_MANAGER_TPP
#define RPC_ANALYZER_MANAGER_TPP

#include "AnalyzerManager.h" 
#include "AnalyzerFactory.h"
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   AnalyzerManager<D>::AnalyzerManager(Simulator<D>& simulator, System<D>& system)
   : Manager< Analyzer<D> >(),
     simulatorPtr_(&simulator),
     systemPtr_(&system)
   {  setClassName("AnalyzerManager"); }

   /*
   * Destructor.
   */
   template <int D>
   AnalyzerManager<D>::~AnalyzerManager()
   {} 
   
   /*
   * Return a pointer to a new AnalyzerFactory object.
   */
   template <int D>
   Factory< Analyzer<D> >* AnalyzerManager<D>::newDefaultFactory() const
   {  return new AnalyzerFactory<D>(*simulatorPtr_, *systemPtr_); }

   /*
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   template <int D>
   void AnalyzerManager<D>::readParameters(std::istream &in)
   {
      Analyzer<D>::baseInterval = 1;
      readOptional(in,"baseInterval", Analyzer<D>::baseInterval);
      Manager< Analyzer<D> >::readParameters(in);
   }

   /*
   * Call initialize method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::setup() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setup();
      }
   }
 
   /*
   * Call sample method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::sample(long iStep) 
   {
      UTIL_CHECK(Analyzer<D>::baseInterval > 0);
      UTIL_CHECK(iStep % Analyzer<D>::baseInterval == 0);
      for (int i = 0; i < size(); ++i) {
         (*this)[i].sample(iStep);
      }
   }
 
   /*
   * Call output method of each analyzer.
   */
   template <int D>
   void AnalyzerManager<D>::output() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].output();
      }
   }

}
}
#endif
