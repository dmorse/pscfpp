#ifndef RPG_ANALYZER_MANAGER_H
#define RPG_ANALYZER_MANAGER_H

#include "Analyzer.h"                  // template parameter
#include <util/param/Manager.h>        // base class template
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class AnalyzerManager : public Manager< Analyzer<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
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
      
      #if 0
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
      
      #endif
      
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
      using ParamComposite::readOptional;
      using ParamComposite::read;
   
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
