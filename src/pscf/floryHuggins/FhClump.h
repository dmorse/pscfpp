#ifndef PSCF_FH_CLUMP_H
#define PSCF_FH_CLUMP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf { 

   using namespace Util;

   /**
   * Sub-unit of a molecule in Flory-Huggins theory.
   *
   * A clump is an aggregate of all monomers of a specified type within 
   * a molecule. A clump has a monomer id and a size. The size of a clump 
   * is the total volume occupied by all monomers of the specified type 
   * in a molecule of a specific species, divided by the monomer reference 
   * volume. 
   * 
   * For a block copolymer, a clump is generally different than a block
   * because a clump may include the monomers in two or more blocks of 
   * the same monomer type. Hompolymer and solvent molecular species 
   * each have only one clump.
   *
   * \ingroup Pscf_FloryHuggins_Module
   */
   class FhClump
   {
   public:
  
      /**
      * Constructor.
      */ 
      FhClump();
  
      /**
      * Serialize to/from archive.
      *
      * \param ar input or output Archive
      * \param versionId archive format version index
      */ 
      template <class Archive>
      void serialize(Archive& ar, unsigned int versionId);

      /// \name Setters
      //@{
    
      /**
      * Set the monomer type id.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size of this clump.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size clump size (volume / monomer volume)
      */ 
      void setSize(double size);
  
      //@}
      /// \name Accessors
      //@{
    
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this block.
      */
      double size() const;

      //@}

   private:
  
      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size = volume / monomer reference volume. 
      double size_;

      friend 
      std::istream& operator >> (std::istream& in, FhClump &block);

      friend 
      std::ostream& operator << (std::ostream& out, const FhClump &block);

   };

   /**
   * istream extractor for a FhClump.
   *
   * \param in  input stream
   * \param block  FhClump to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, FhClump &block);

   /**
   * ostream inserter for a FhClump.
   *
   * \param out  output stream
   * \param block  FhClump to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const FhClump &block);

   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   inline int FhClump::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double FhClump::size() const
   {  return size_; }
    
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void FhClump::serialize(Archive& ar, unsigned int)
   {
      ar & monomerId_;
      ar & size_;
   }
    
} 
#endif 
