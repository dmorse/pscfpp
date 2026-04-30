#ifndef RPC_COMPRESSOR_H
#define RPC_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>

namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class Compressor : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Compressor();

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      Compressor(System<D>& system);

      /**
      * Destructor.
      */
      virtual ~Compressor() = default;

      /**
      * Iterate Langrange multiplier field.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int compress() = 0;

      /**
      * Output report of timing results to stream.
      *
      * \param out output stream for results
      */
      virtual void outputTimers(std::ostream& out) const = 0;

      /**
      * Clear timers.
      */
      virtual void clearTimers() = 0;

      /**
      * Get the number of times the MDE has been solved.
      */
      int mdeCounter() const;

   protected:

      /**
      *  Create association with the parent System.
      *
      * \param system parent System object
      */
      void setSystem(System<D>& system);

      /**
      * Return parent system by const reference.
      */
      System<D> const & system() const;

      /**
      * Return parent system by non-const reference.
      */
      System<D>& system();

      /**
      * Count how many times MDE has been solved.
      */
      int mdeCounter_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Inline member functions
 
   /*
   * Return parent system by const reference.
   */
   template <int D> inline
   System<D> const & Compressor<D>::system() const
   {  
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

   /*
   * Return parent system by non-const reference.
   */
   template <int D> inline
   System<D>& Compressor<D>::system()
   {  
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

   // Non-inline Member functions

   /*
   * Default constructor.
   */
   template <int D>
   Compressor<D>::Compressor()
    : mdeCounter_(0),
      sysPtr_(nullptr)
   {  setClassName("Compressor"); }

   /*
   * Constructor (creates association with parent system)
   */
   template <int D>
   Compressor<D>::Compressor(System<D>& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   /*
   * Create association with the parent system.
   */
   template <int D>
   void Compressor<D>::setSystem(System<D>& system)
   {  sysPtr_ = &system; }

   /*
   * Get number of times MDE has been solved.
   */
   template <int D>
   inline int Compressor<D>::mdeCounter() const
   {  return mdeCounter_; }

} // namespace Rpc
} // namespace Pscf
#endif
