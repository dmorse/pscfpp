#ifndef UTIL_BIT_H
#define UTIL_BIT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Represents a specific bit location within an unsigned int.
   *
   * Provides methods to query, set or clear a particular bit.
   *
   * \ingroup Misc_Module
   */
   class Bit 
   {
   public:
  
      /**
      * Default constructor.
      */ 
      Bit();
  
      /**
      * Constructor.
      *
      * \param shift location of the bit, 0 < shift <= 32.
      */ 
      Bit(unsigned int shift);
  
      /**
      * Set or reset the bit mask.
      *
      * \param shift location of the bit, 0 < shift <= 32.
      */ 
      void setMask(unsigned int shift);
  
      /**
      * Set this bit in the flags parameter
      *
      * \param flags unsigned int to be modified
      */
      void set(unsigned int& flags) const;
   
      /**
      * Clear this bit in the flags parameter
      *
      * \param flags unsigned int to be modified
      */
      void clear(unsigned int& flags) const;

      /**
      * Is this bit set in the flags integer?
      *
      * \param flags unsigned int to be queried
      */
      bool isSet(unsigned int flags) const;
   
      /**
      * Return integer with only this bit set.
      */
      unsigned int mask() const;

   private:

      /// Integer with this bit set, all others clear.
      unsigned int mask_;
   
   };

   /*
   * Set this bit in the flags parameter.
   */
   inline void Bit::set(unsigned int& flags) const
   {  flags |= mask_; }

   /*
   * Clear this bit in the flags parameter.
   */
   inline void Bit::clear(unsigned int& flags) const
   {  flags &= (~mask_); }

   /*
   * Is this bit set in the flags integer?
   */
   inline bool Bit::isSet(unsigned int flags) const
   {  return bool(flags & mask_); }

   /*
   * Return unsigned int with only this bit set.
   */
   inline unsigned int Bit::mask() const
   {  return mask_; }
   
}
#endif
