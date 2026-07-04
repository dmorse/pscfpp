#ifndef RPG_ANALYZER_FACTORY_H
#define RPG_ANALYZER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/analyzer/Analyzer.h>
//#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Factory for subclasses of Analyzer.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class AnalyzerFactory : public Factory< Rp::Analyzer<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rp::Simulator<D, Rpg::Types<D> > object
      * \param system  parent Rp::System<D, Rpg::Types<D> > object
      */
      AnalyzerFactory(Rp::Simulator<D, Rpg::Types<D> >& simulator, Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Method to create any Analyzer supplied with PSCF.
      *
      * \param className name of the Analyzer subclass
      * \return Analyzer* pointer to new instance of className
      */
      Rp::Analyzer<D, Rpg::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::Analyzer<D, Rpg::Types<D> > >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Rp::Simulator<D, Rpg::Types<D> >* simPtr_;

      /// Pointer to the parent system.
      Rp::System<D, Rpg::Types<D> >* sysPtr_;
      
   };

   // Explicit instantiation declarations
   extern template class AnalyzerFactory<1>;
   extern template class AnalyzerFactory<2>;
   extern template class AnalyzerFactory<3>;

}
}
#endif
