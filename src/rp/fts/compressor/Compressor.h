#ifndef RP_COMPRESSOR_H
#define RP_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {


   using namespace Util;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T>
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
      Compressor(System<D,T>& system);

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
      void setSystem(System<D,T>& system);

      /**
      * Return parent system by const reference.
      */
      System<D,T> const & system() const;

      /**
      * Return parent system by non-const reference.
      */
      System<D,T>& system();

      /**
      * Count how many times MDE has been solved.
      */
      int mdeCounter_;

   private:

      /// Pointer to the associated system object.
      System<D,T>* sysPtr_;

   };

   // Member functions (all inline)
 
   /*
   * Default constructor.
   */
   template <int D, class T> inline
   Compressor<D,T>::Compressor()
    : mdeCounter_(0),
      sysPtr_(nullptr)
   {  setClassName("Compressor"); }

   /*
   * Constructor (creates association with parent system)
   */
   template <int D, class T> inline
   Compressor<D,T>::Compressor(System<D,T>& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   /*
   * Create association with the parent system.
   */
   template <int D, class T> inline
   void Compressor<D,T>::setSystem(System<D,T>& system)
   {  sysPtr_ = &system; }

   /*
   * Get number of times MDE has been solved.
   */
   template <int D, class T> inline
   int Compressor<D,T>::mdeCounter() const
   {  return mdeCounter_; }

   /*
   * Return parent system by const reference.
   */
   template <int D, class T> inline
   System<D,T> const & Compressor<D,T>::system() const
   {  
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

   /*
   * Return parent system by non-const reference.
   */
   template <int D, class T> inline
   System<D,T>& Compressor<D,T>::system()
   {  
      UTIL_ASSERT(sysPtr_);
      return *sysPtr_;
   }

} // namespace Rp
} // namespace Pscf
#endif
