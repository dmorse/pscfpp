#ifndef RPC_ANALYZER_H
#define RPC_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>  // base class template
#include <rpc/system/Types.h>          // base class template argument

namespace Pscf {
namespace Rpc {


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
   * their entire public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::Analyzer
   * \see \ref psfts_algo_analyzer_page "Analyzers"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class Analyzer : public Rp::Analyzer<D, Rp::Simulator<D, Rpc::Types<D> >, Rp::System<D, Rpc::Types<D> > >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      Analyzer(Rp::Simulator<D, Rpc::Types<D> >& simulator, Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Destructor.
      */
      virtual ~Analyzer() = default;

      Analyzer() = delete;
      Analyzer(Analyzer<D> const &) = delete;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template 
      class Analyzer<1, Rp::Simulator<1, Rpc::Types<1> >, Rp::System<1, Rpc::Types<1> > >;
      extern template 
      class Analyzer<2, Rp::Simulator<2, Rpc::Types<2> >, Rp::System<2, Rpc::Types<2> > >;
      extern template 
      class Analyzer<3, Rp::Simulator<3, Rpc::Types<3> >, Rp::System<3, Rpc::Types<3> > >;
   } 
   namespace Rpc {
      extern template class Analyzer<1>;
      extern template class Analyzer<2>;
      extern template class Analyzer<3>;
   } 
} 
#endif
