#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include "Analyzer.h"                  // base class template parameter
#include <util/param/Manager.h>        // base class template

namespace Pscf {
namespace Rpc {

   using namespace Util;
   
   template <int D> class System;
   template <int D> class Simulator;


   /**
   * Manager for a list of Analyzer objects.
   *
   * \see rpc_AnalyzerManager_page Manual Page
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AnalyzerManager : public Manager< Analyzer<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Simulator
      * \param system parent System
      */
      AnalyzerManager(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~AnalyzerManager();

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      virtual void readParameters(std::istream &in);
      
      /**
      * Call initialize method of each Analyzer.
      * 
      * This method should be called just before the main
      * simulation loop, after an initial configuration is
      * known. It calls the setup() functionfor each 
      * analyzer, or does nothing if size() == 0.
      */
      void setup();
 
      /**
      * Call sample method of each Analyzer.
      *
      * \pre Analyzer::baseInterval > 0
      * \pre iStep::baseInterval == 0
      * 
      * \param iStep step counter for main loop
      */
      void sample(long iStep);
 
      /**
      * Call output method of each analyzer.
      * 
      * This method should be called after the main
      * simulation loop. It calls the output() function
      * of each analyzer, or does nothing if size() == 0.
      */
      void output();

      using Manager< Analyzer<D> >::size;

   protected:

      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
   
   private:

      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;
      
      /**
      * Pointer to parent System.
      */
      System<D>* systemPtr_;

      /**
      * Return pointer to a new AnalyzerFactory.
      */
      virtual Factory< Analyzer<D> >* newDefaultFactory() const;

   };

   // Explicit instantiation declarations
   extern template class AnalyzerManager<1>;
   extern template class AnalyzerManager<2>;
   extern template class AnalyzerManager<3>;

}
}
#endif
