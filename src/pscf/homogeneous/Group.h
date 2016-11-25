#ifndef PSCF_HOMOGENEOUS_GROUP_H
#define PSCF_HOMOGENEOUS_GROUP_H

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
   * A linear homopolymer block within a block copolymer.
   *
   * \ingroup Pscf_Homogeneous_Module
   */
   class Group
   {
   public:
  
      /**
      * Constructor.
      */ 
      Group();
  
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
      std::istream& operator >> (std::istream& in, Group &block);

      friend 
      std::ostream& operator << (std::ostream& out, const Group &block);

   };

   /**
   * istream extractor for a Group.
   *
   * \param in  input stream
   * \param block  Group to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Group &block);

   /**
   * ostream inserter for a Group.
   *
   * \param out  output stream
   * \param block  Group to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Group &block);

   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   inline int Group::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double Group::size() const
   {  return size_; }
    
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Group::serialize(Archive& ar, unsigned int)
   {
      ar & monomerId_;
      ar & size_;
   }
    
}
} 
#endif 
