#ifndef RP_ANALYZER_MANAGER_H
#define RP_ANALYZER_MANAGER_H

#include <util/param/Manager.h>        // base class template

// Forward declaration
namespace Util {
   template <class T> class Factory;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named AnalyzerManager,
   * that are defined in Rpc and Rpg namespaces.
   *
   * Template parameters:
   *
   *   - D : dimension of space (1, 2, or 3)
   *   - T : Types class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \see \ref rp_AnalyzerManager_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class AnalyzerManager : public Manager< typename T::Analyzer >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator
      * \param system  parent System
      */
      AnalyzerManager(typename T::Simulator& simulator, 
                      typename T::System& system);

      /**
      * Destructor.
      */
      virtual ~AnalyzerManager();

      /**
      * Read body of parameter file block.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Call the setup function of each Analyzer.
      *
      * This function should be called just before the main
      * simulation loop, after an initial configuration is
      * known. It calls the setup() functionfor each
      * analyzer, or does nothing if size() == 0.
      */
      void setup();

      /**
      * Call the sample function of each Analyzer.
      *
      * \pre T::Analyzer::baseInterval > 0
      * \pre iStep % T::Analyzer::baseInterval == 0
      *
      * \param iStep  step counter for main loop
      */
      void sample(long iStep);

      /**
      * Call the output function of each analyzer.
      *
      * This function should be called after the main
      * simulation loop. It calls the output() function
      * of each analyzer, or does nothing if size() == 0.
      */
      void output();

   private:

      using AnalyzerT = typename T::Analyzer;
      using AnalyzerFactoryT = typename T::AnalyzerFactory;
      using Base = Manager< AnalyzerT >;

      /**
      * Pointer to parent Simulator
      */
      typename T::Simulator* simulatorPtr_;

      /**
      * Pointer to parent System.
      */
      typename T::System* systemPtr_;

      /**
      * Return pointer to a new AnalyzerFactory.
      */
      Factory<typename T::Analyzer>* newDefaultFactory() const override;

   };

}
}
#endif
