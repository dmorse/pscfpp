#ifndef RPG_ANALYZER_H
#define RPG_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>   // class template
#include <rpg/system/Types.h>           // class argument

#if 0
namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * The periodic action associated with an Analyzer may involve retrieval
   * or computation of a physical property value, adding it to a statistical
   * accumulator, and/or outputting it to file. This periodic action must
   * be implemented by the pure virtual sample() function.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::Analyzer, and inherit
   * their public entire interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::Analyzer
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class Analyzer : public Rp::Analyzer<D, Rpg::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      Analyzer(Rp::Simulator<D, Rpg::Types<D> >& simulator, 
               Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Analyzer() = default;

      using Rp::Analyzer<D, Rpg::Types<D> >::baseInterval;

   };

} // namespace Rpg
} // namespace Pscf
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, Rpg::Types<1> >;
      extern template class Analyzer<2, Rpg::Types<2> >;
      extern template class Analyzer<3, Rpg::Types<3> >;
   } 
} 
#endif
