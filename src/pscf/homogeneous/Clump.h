#ifndef PSCF_HOMOGENEOUS_CLUMP_H
#define PSCF_HOMOGENEOUS_CLUMP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf { 
namespace Homogeneous { 

   using namespace Util;

   /**
   * Set of all monomers of the same type in a molecule.
   *
   * A clump has a monomer id and a size. The size of a clump is 
   * the volume occupied by all monomers of the specified type in
   * particular molecular species, divided by a monomer reference 
   * volume. 
   * 
   * For a block copolymer, a clump is generally different than 
   * a block because a clump may include the monomers in two or 
   * more blocks of the same monomer type. Hompolymer and point
   * solvent molecular species each have only one clump.
   *
   * \ingroup Pscf_Homogeneous_Module
   */
   class Clump
   {
   public:
  
      /**
      * Constructor.
      */ 
      Clump();
  
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
      * Set the monomer id.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size of this block.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size block size (number of monomers).
      */ 
      void setSize(double size);
  
      //@}
      /// \name Accessors (getters)
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
      std::istream& operator >> (std::istream& in, Clump &block);

      friend 
      std::ostream& operator << (std::ostream& out, const Clump &block);

   };

   /**
   * istream extractor for a Clump.
   *
   * \param in  input stream
   * \param block  Clump to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Clump &block);

   /**
   * ostream inserter for a Clump.
   *
   * \param out  output stream
   * \param block  Clump to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Clump &block);

   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   inline int Clump::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double Clump::size() const
   {  return size_; }
    
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Clump::serialize(Archive& ar, unsigned int)
   {
      ar & monomerId_;
      ar & size_;
   }
    
}
} 
#endif 
