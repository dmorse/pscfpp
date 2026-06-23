#ifndef RPG_COMPRESSOR_H
#define RPG_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/system/Types.h>           // template argument
#include <util/param/ParamComposite.h>  // base class
#include <util/global.h>

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Rpg_Fts_Compressor_Module
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
      Compressor(Rp::System<D, Rpg::Types<D> >& system);

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
      * Get the number of times the MDE has been solved.
      */
      int mdeCounter() const;

      /**
      * Log output timing results.
      *
      * \param out  output stream
      */
      virtual void outputTimers(std::ostream& out) const = 0;

      /**
      * Clear timers
      */
      virtual void clearTimers() = 0;

   protected:

      /**
      * Count how many times MDE has been solved.
      */
      int mdeCounter_;

      /**
      *  Create association with the parent System.
      *
      * \param system  parent System object
      */
      void setSystem(Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Return parent system by const reference.
      */
      Rp::System<D, Rpg::Types<D> > const & system() const;

      /**
      * Return parent system by non-const reference.
      */
      Rp::System<D, Rpg::Types<D> >& system();

   private:

      /// Pointer to the associated system object.
      Rp::System<D, Rpg::Types<D> >* sysPtr_;

   };

   // Inline member functions

   /*
   * Return parent system by const reference.
   */
   template <int D> inline
   Rp::System<D, Rpg::Types<D> > const & Compressor<D>::system() const
   {
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

   /*
   * Return parent system by non-const reference.
   */
   template <int D> inline
   Rp::System<D, Rpg::Types<D> >& Compressor<D>::system()
   {
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

   // Non-inline member functions

   /*
   * Default constructor.
   */
   template <int D>
   Compressor<D>::Compressor()
    : mdeCounter_(0),
      sysPtr_(nullptr)
   {  setClassName("Compressor"); }

   /*
   * Constructor.
   */
   template <int D>
   Compressor<D>::Compressor(Rp::System<D, Rpg::Types<D> >& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   /*
   * Get number of times MDE has been solved.
   */
   template <int D>
   int Compressor<D>::mdeCounter() const
   {  return mdeCounter_; }

   // Protected function

   /*
   * Create association with the parent system.
   */
   template <int D>
   void Compressor<D>::setSystem(Rp::System<D, Rpg::Types<D> >& system)
   {
      UTIL_CHECK(!sysPtr_);
      sysPtr_ = &system;
   }

} // namespace Rpg
} // namespace Pscf
#endif
