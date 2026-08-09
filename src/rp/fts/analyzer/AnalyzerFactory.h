#ifndef RP_ANALYZER_FACTORY_H
#define RP_ANALYZER_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>         // base class template
#include <rp/fts/analyzer/Analyzer.h>   // base class argument

#include <pscf/backends/TmplDeclare.h>
#include <string>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Factory for subclasses of Analyzer.
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class AnalyzerFactory : public Factory< Analyzer<D,T> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AnalyzerFactory(Simulator<D,T>& simulator, 
		      System<D,T>& system);

      /**
      * Create an instance of a specified subclass of Analyzer..
      *
      * \param className  name of the desired subclass of Analyzer<D,T>
      * \return pointer to new instance of the specified class
      */
      Analyzer<D,T>* factory(const std::string &className) const;

      using Factory< Analyzer<D,T> >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Simulator<D,T>* simPtr_;

      /// Pointer to the parent system.
      System<D,T>* sysPtr_;
      
   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(AnalyzerFactory)

} // namespace Rp
} // namespace Pscf
#endif
