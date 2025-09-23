#ifndef RPG_COMPRESSOR_H
#define RPG_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Rpg {

   template <int D> class System;

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
      Compressor(System<D>& system);

      /**
      * Destructor.
      */
      ~Compressor();

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
      * Log output timing results 
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
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  return *sysPtr_;}
      
      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}
      
   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   #if 0
   // Default constructor
   template <int D>
   inline Compressor<D>::Compressor()
    : sysPtr_(&system)
   {  setClassName("Compressor"); }
   #endif

   // Constructor
   template <int D>
   Compressor<D>::Compressor(System<D>& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   // Destructor
   template <int D>
   Compressor<D>::~Compressor()
   {}

   // Get number of times MDE has been solved.
   template <int D>
   inline int Compressor<D>::mdeCounter() const
   {  return mdeCounter_; }

} // namespace Rpg
} // namespace Pscf
#endif
